import os
import sys
import argparse
import shutil
import logging
from datetime import datetime
from typing import List, Dict

# Add src directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.config import ConfigManager
from src.database import DatabaseInterface
from src.search_strategies import BindingSiteSearcher
from src.filtering import SequenceFilter
from src.utils import load_gene_list


def setup_logging_and_copy_config(config_file: str, output_dir: str):
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
    
    # Copy config file to configs subdirectory
    config_dest = os.path.join(configs_dir, "config_used.yaml")
    shutil.copy2(config_file, config_dest)
    logging.info(f"Configuration copied to: {config_dest}")
    
    # Copy species config if it exists
    species_config_path = os.path.join(os.path.dirname(config_file), "species_config.json")
    if os.path.exists(species_config_path):
        species_dest = os.path.join(configs_dir, "species_config.json")
        shutil.copy2(species_config_path, species_dest)
        logging.info(f"Species configuration copied to: {species_dest}")
    
    return log_file




def main():
    parser = argparse.ArgumentParser(description="Find isoform-specific binding sites per gene using Ensembl and IsoformSpecificStrategy")
    parser.add_argument('--genes_file', required=True, help='Path to a text file of gene names (one per line)')
    parser.add_argument('--config', default='configs/config_specific.yaml', help='Path to configuration file')
    parser.add_argument('--species', default='mouse', help='Species name (mouse, human, rat, zebrafish)')
    parser.add_argument('--output_dir', help='Override output directory from config')
    parser.add_argument('--genome_fasta', help='Override genome FASTA path from config')
    parser.add_argument('--blast_species', nargs='+', help='Target species for BLAST search (default: from config blast.species)')
    args = parser.parse_args()

    # Load configuration
    cfg = ConfigManager(args.config, species=args.species)
    
    # Override config with command line arguments if provided
    if args.output_dir:
        cfg.output.output_dir = args.output_dir
    if args.genome_fasta:
        cfg.genome.genome_fasta_path = args.genome_fasta
    if args.blast_species:
        cfg.blast.species = args.blast_species

    # Validate configuration
    errors = cfg.validate_config()
    if errors:
        print("Configuration errors:")
        for error in errors:
            print(f"  - {error}")
        return
    
    # Setup logging and copy configuration files
    setup_logging_and_copy_config(args.config, cfg.output.output_dir)

    # Initialize database interface
    db = DatabaseInterface(cfg.database)

    # Prepare genome accessor - prefer local, fallback to Ensembl
    accessor = None
    if cfg.genome.use_local_first and cfg.genome.genome_fasta_path and os.path.exists(cfg.genome.genome_fasta_path):
        print(f"Using local genome: {cfg.genome.genome_fasta_path}")
        accessor = db.local_genome_accessor(cfg.genome.genome_fasta_path)
    
    if accessor is None:
        print("Using Ensembl genome accessor")
        accessor = db.ensembl_genome_accessor()

    # Read gene list
    genes = load_gene_list(args.genes_file)
    print(f"Processing {len(genes)} genes...")

    # Build sequences and isoforms maps
    sequences: Dict[str, Dict] = {}
    isoforms_map: Dict[str, List[Dict]] = {}
    for g in genes:
        isoforms = db.get_normalized_isoforms(g)
        if not isoforms:
            print(f"Warning: No isoforms found for gene {g}")
            continue
        isoforms_map[g] = isoforms
        sequences[g] = {'sequence': isoforms[0].get('seq', 'N' * cfg.search.binding_site_length)}

    print(f"Found isoforms for {len(isoforms_map)} genes")

    # Initialize searcher
    searcher = BindingSiteSearcher(cfg.search, cfg.filter, genome_accessor=accessor)
    binding_sites = searcher.search_all_genes(sequences, isoforms=isoforms_map)

    # Create output directory first
    os.makedirs(cfg.output.output_dir, exist_ok=True)

    # Initialize filter for BLAST analysis
    filter = SequenceFilter(cfg.filter, cfg.blast)

    print("Starting BLAST analysis...")

    # 1. Pre-BLAST filtering
    print("Applying pre-BLAST filters...")
    filtered_sites = filter.pre_blast_filter(binding_sites)

    # 2. Run BLAST analysis
    print("Running BLAST analysis...")
    blast_results_file = filter.run_blast(filtered_sites, cfg.output.output_dir)

    # 3. Specificity filtering
    print("Applying specificity filters...")
    final_sites = filter.specificity_filter(filtered_sites, blast_results_file)

    # 4. Save final results
    print("Saving final results...")
    filter.save_filtered_results(final_sites, cfg.output.output_dir)

    # Also save original results for comparison
    os.makedirs(cfg.output.output_dir, exist_ok=True)
    searcher.save_binding_sites(binding_sites, cfg.output.output_dir)
    
    print(f"Results saved to {cfg.output.output_dir}")
    print(f"Original sites: {sum(len(sites) for sites in binding_sites.values())}")
    print(f"After pre-BLAST filter: {sum(len(sites) for sites in filtered_sites.values())}")
    print(f"After specificity filter: {sum(len(sites) for sites in final_sites.values())}")


if __name__ == '__main__':
    main()
