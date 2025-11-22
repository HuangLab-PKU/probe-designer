import os
import sys
import argparse
import json
import shutil
import logging
from datetime import datetime
from typing import List, Dict

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.config import ConfigManager
from src.database import DatabaseInterface
from src.search_strategies import BindingSiteSearcher
from src.filtering import SequenceFilter
from src.utils import load_gene_list


def setup_logging_and_copy_config(config_file: str, output_dir: str, cfg=None):
    """Setup logging and copy configuration files to output directory."""
    # Create logs directory in output
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    
    # Create configs directory in output
    configs_dir = os.path.join(output_dir, "configs")
    os.makedirs(configs_dir, exist_ok=True)
    
    # Setup logging
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(logs_dir, f"probe_design_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    # Save the actual configuration used (with all overrides) to configs subdirectory
    config_dest = os.path.join(configs_dir, "config_used.yaml")
    if cfg is not None:
        # Save the actual configuration object with all overrides
        cfg.save_config(config_dest)
        logging.info(f"Actual configuration (with overrides) saved to: {config_dest}")
    else:
        # Fallback: copy original config file
        shutil.copy2(config_file, config_dest)
        logging.info(f"Original configuration copied to: {config_dest}")
    
    # Copy species config if it exists
    species_config_path = os.path.join(os.path.dirname(config_file), "species_config.json")
    if os.path.exists(species_config_path):
        species_dest = os.path.join(configs_dir, "species_config.json")
        shutil.copy2(species_config_path, species_dest)
        logging.info(f"Species configuration copied to: {species_dest}")
    
    return log_file


def main():
    parser = argparse.ArgumentParser(description="Find consensus binding sites (max isoforms per gene) using Ensembl and IsoformConsensusStrategy")
    parser.add_argument('--config', default='configs/config_consensus.yaml', help='Path to configuration file')
    parser.add_argument('--genes_file', help='Path to a text file of gene names (one per line, overrides config)')
    parser.add_argument('--species', help='Override species name (mouse, human, rat, zebrafish)')
    parser.add_argument('--blast_species', nargs='+', help='Override target species for BLAST search (default: from config blast.species)')
    parser.add_argument('--genome_fasta', help='Override genome FASTA path from config')
    parser.add_argument('--output_dir', help='Override output directory from config')
    parser.add_argument('--progress', action='store_true', help='Show progress bar during binding site search')
    args = parser.parse_args()

    # Load configuration
    cfg = ConfigManager(args.config)
    
    # Override config with command line arguments if provided
    if args.species:
        cfg.species_config.species = args.species
    if args.blast_species:
        cfg.blast.species = args.blast_species
    if args.genome_fasta:
        cfg.genome.genome_fasta_path = args.genome_fasta
    if args.output_dir:
        cfg.output.output_dir = args.output_dir
        print(f"Overriding output_dir with command line argument: {args.output_dir}")
    else:
        print(f"Using output_dir from config: {cfg.output.output_dir}")
    
    # Get genes_file from command line or config
    genes_file = args.genes_file or cfg.search.genes_file
    if not genes_file:
        parser.error("genes_file must be specified either in config file (search.genes_file) or via --genes_file argument")
    
    if args.genes_file:
        print(f"Using genes_file from command line: {args.genes_file}")
    else:
        print(f"Using genes_file from config: {genes_file}")

    # Validate configuration
    errors = cfg.validate_config()
    if errors:
        print("Configuration errors:")
        for error in errors:
            print(f"  - {error}")
        return
    
    # Setup logging and copy configuration files (after all overrides)
    setup_logging_and_copy_config(args.config, cfg.output.output_dir, cfg)

    # Initialize database interface
    db = DatabaseInterface(cfg.database)

    # Initialize genome accessor
    accessor = db.initialize_genome_accessor(cfg.genome)

    # Read gene list
    genes = load_gene_list(genes_file)
    print(f"Processing {len(genes)} genes...")

    # Create output directory
    os.makedirs(cfg.output.output_dir, exist_ok=True)
    
    # Create cache directories for different types of intermediate results
    print("Building gene-isoforms map...")
    isoform_info_dir = os.path.join(cfg.output.output_dir, "isoform_info")
    os.makedirs(isoform_info_dir, exist_ok=True)    

    # build gene-isoforms map
    isoforms_map: Dict[str, List[Dict]] = {}
    cached_genes = 0
    
    for g in genes:
        isoform_file = os.path.join(isoform_info_dir, f"{g}.json")
        if os.path.exists(isoform_file):
            # Load cached isoform data for this gene
            with open(isoform_file, 'r', encoding='utf-8') as f:
                isoforms = json.load(f)
            cached_genes += 1
        else:
            # build isoforms for this gene
            isoforms = db.get_isoform_info(g)
            # save to cache (one file per gene)
            with open(isoform_file, 'w', encoding='utf-8') as f:
                json.dump(isoforms, f, indent=2, ensure_ascii=False)      
        isoforms_map[g] = isoforms
    print(f"Found isoforms for {len(isoforms_map)} genes ({cached_genes} from cache)")

    # search binding sites for each gene
    print("Searching binding sites...")
    binding_site_candidates_dir = os.path.join(cfg.output.output_dir, "binding_site_candidates")
    os.makedirs(binding_site_candidates_dir, exist_ok=True)

    binding_sites: Dict[str, List[Dict]] = {}
    cached_binding_genes = 0

    # initialize searcher
    searcher = BindingSiteSearcher(cfg.search, cfg.filter, genome_accessor=accessor, show_progress=args.progress)
    
    for g in genes:
        binding_site_file = os.path.join(binding_site_candidates_dir, f"{g}.json")
        if os.path.exists(binding_site_file):
            # Load cached binding sites for this gene
            with open(binding_site_file, 'r', encoding='utf-8') as f:
                gene_sites = json.load(f)
            cached_binding_genes += 1
        else:
            # Search for binding sites for this gene
            if g not in isoforms_map:
                print(f"Warning: No isoforms found for gene {g}")
                continue
            gene_isoforms = {g: isoforms_map[g]}
            gene_binding_sites = searcher.search_all_genes(sequences=None, isoforms=gene_isoforms)
            gene_sites = gene_binding_sites[g]
            # save to cache (one file per gene)
            with open(binding_site_file, 'w', encoding='utf-8') as f:
                json.dump(gene_sites, f, indent=2, ensure_ascii=False)
        
        if gene_sites:
            binding_sites[g] = gene_sites

    print(f"Found binding sites for {len(binding_sites)} genes ({cached_binding_genes} from cache)")
    
    # Save combined binding sites for compatibility
    searcher.save_binding_sites(binding_sites, cfg.output.output_dir)
    
    # Initialize filter for BLAST analysis
    filter = SequenceFilter(cfg.filter, cfg.blast)
    
    print("Starting BLAST analysis...")
    blast_results_dir = os.path.join(cfg.output.output_dir, "blast_results")
    os.makedirs(blast_results_dir, exist_ok=True)

    # 1. Run BLAST analysis
    print("Running BLAST analysis...")
    blast_results_file = filter.run_blast(binding_sites, blast_results_dir)
    
    # 2. Specificity filtering
    print("Applying specificity filters...")
    final_sites = filter.specificity_filter(binding_sites, blast_results_file)
    
    # 3. Save final results
    print("Saving final results...")
    filter.save_filtered_results(final_sites, cfg.output.output_dir)
    
    print(f"Results saved to {cfg.output.output_dir}")
    print(f"Original sites: {sum(len(sites) for sites in binding_sites.values())}")
    print(f"After specificity filter: {sum(len(sites) for sites in final_sites.values())}")


if __name__ == '__main__':
    main()
