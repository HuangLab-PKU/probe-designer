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


def read_gene_list(path: str) -> List[str]:
    """Read gene list from text file or Excel file."""
    if path.endswith('.xlsx') or path.endswith('.xls'):
        import pandas as pd
        df = pd.read_excel(path)
        # Assume first column contains gene names
        return [str(gene).strip() for gene in df.iloc[:, 0] if str(gene).strip() and str(gene).strip() != 'nan']
    else:
        with open(path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip()]


def setup_logging_and_copy_config(config_file: str, output_dir: str, cfg=None):
    """Setup logging and copy configuration files to output directory."""
    # Create logs directory in output
    logs_dir = os.path.join(output_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    
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
    
    # Save the actual configuration used (with all overrides) to output directory
    config_dest = os.path.join(output_dir, "config_used.yaml")
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
        species_dest = os.path.join(output_dir, "species_config.json")
        shutil.copy2(species_config_path, species_dest)
        logging.info(f"Species configuration copied to: {species_dest}")
    
    return log_file


def build_isoforms_map(db: DatabaseInterface, gene: str) -> List[Dict]:
    # Use Ensembl isoform details
    txs = db.get_isoform_sequences(gene)
    isoforms = []
    for tx in txs:
        exons = tx.get('exons', [])
        if not exons:
            continue
        starts = [e.get('start') for e in exons if 'start' in e]
        ends = [e.get('end') for e in exons if 'end' in e]
        if not starts or not ends:
            continue
        isoforms.append({
            'id': tx.get('id'),
            'external_name': tx.get('external_name') or tx.get('id'),
            'exons': [{'start': e['start'], 'end': e['end']} for e in exons if 'start' in e and 'end' in e],
            'start': min(starts),
            'end': max(ends),
            'strand': tx.get('strand', 1),
            'seq_region_name': tx.get('seq_region_name', '1'),
            'seq': tx.get('seq', '')
        })
    return isoforms


def main():
    parser = argparse.ArgumentParser(description="Find consensus binding sites (max isoforms per gene) using Ensembl and IsoformConsensusStrategy")
    parser.add_argument('--genes_file', required=True, help='Path to a text file of gene names (one per line)')
    parser.add_argument('--config', default='configs/config_consensus.yaml', help='Path to configuration file')
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
    
    # Setup logging and copy configuration files (after all overrides)
    setup_logging_and_copy_config(args.config, cfg.output.output_dir, cfg)

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
    genes = read_gene_list(args.genes_file)
    print(f"Processing {len(genes)} genes...")

    # Create output directory first
    os.makedirs(cfg.output.output_dir, exist_ok=True)
    
    # Create cache directories for different types of intermediate results
    isoform_info_dir = os.path.join(cfg.output.output_dir, "isoform_info")
    binding_site_candidates_dir = os.path.join(cfg.output.output_dir, "binding_site_candidates")
    blast_results_dir = os.path.join(cfg.output.output_dir, "blast_results")
    
    os.makedirs(isoform_info_dir, exist_ok=True)
    os.makedirs(binding_site_candidates_dir, exist_ok=True)
    os.makedirs(blast_results_dir, exist_ok=True)
    
    # Check for cached isoforms data (one file per gene)
    print("Loading or building isoforms data...")
    sequences: Dict[str, Dict] = {}
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
            # Build isoforms map for this gene
            isoforms = build_isoforms_map(db, g)
            if not isoforms:
                print(f"Warning: No isoforms found for gene {g}")
                continue
            
            # Save to cache (one file per gene)
            with open(isoform_file, 'w', encoding='utf-8') as f:
                json.dump(isoforms, f, indent=2, ensure_ascii=False)
        
        isoforms_map[g] = isoforms
        # Provide a placeholder sequence (first isoform seq) to satisfy searcher input
        sequences[g] = {'sequence': isoforms[0].get('seq', 'N' * cfg.search.binding_site_length)}

    print(f"Found isoforms for {len(isoforms_map)} genes ({cached_genes} from cache)")

    # Check for cached binding sites (one file per gene)
    print("Loading or searching for binding sites...")
    binding_sites: Dict[str, List[Dict]] = {}
    cached_binding_genes = 0
    
    for g in genes:
        binding_site_file = os.path.join(binding_site_candidates_dir, f"{g}.json")
        
        if os.path.exists(binding_site_file):
            # Load cached binding sites for this gene
            with open(binding_site_file, 'r', encoding='utf-8') as f:
                gene_sites = json.load(f)
            cached_binding_genes += 1
        else:
            # Search for binding sites for this gene
            if g not in sequences:
                continue
            searcher = BindingSiteSearcher(cfg.search, cfg.filter, genome_accessor=accessor)
            gene_sequences = {g: sequences[g]}
            gene_isoforms = {g: isoforms_map[g]}
            gene_binding_sites = searcher.search_all_genes(gene_sequences, isoforms=gene_isoforms)
            gene_sites = gene_binding_sites.get(g, [])
            
            # Save to cache (one file per gene)
            with open(binding_site_file, 'w', encoding='utf-8') as f:
                json.dump(gene_sites, f, indent=2, ensure_ascii=False)
        
        if gene_sites:  # Only add if we have sites
            binding_sites[g] = gene_sites

    print(f"Found binding sites for {len(binding_sites)} genes ({cached_binding_genes} from cache)")
    
    # Save combined binding sites for compatibility
    searcher = BindingSiteSearcher(cfg.search, cfg.filter, genome_accessor=accessor)
    searcher.save_binding_sites(binding_sites, cfg.output.output_dir)
    
    # Initialize filter for BLAST analysis
    filter = SequenceFilter(cfg.filter, cfg.blast)
    
    print("Starting BLAST analysis...")
    
    # 1. Pre-BLAST filtering
    print("Applying pre-BLAST filters...")
    filtered_sites = filter.pre_blast_filter(binding_sites)
    
    # Save pre-BLAST filtered results
    pre_blast_cache_file = os.path.join(cfg.output.output_dir, "pre_blast_filtered.json")
    with open(pre_blast_cache_file, 'w', encoding='utf-8') as f:
        json.dump(filtered_sites, f, indent=2, ensure_ascii=False)
    print("Saved pre-BLAST filtered results")
    
    # 2. Run BLAST analysis
    print("Running BLAST analysis...")
    blast_results_file = filter.run_blast(filtered_sites, blast_results_dir)
    
    # 3. Specificity filtering
    print("Applying specificity filters...")
    final_sites = filter.specificity_filter(filtered_sites, blast_results_file)
    
    # 4. Save final results
    print("Saving final results...")
    filter.save_filtered_results(final_sites, cfg.output.output_dir)
    
    print(f"Results saved to {cfg.output.output_dir}")
    print(f"Original sites: {sum(len(sites) for sites in binding_sites.values())}")
    print(f"After pre-BLAST filter: {sum(len(sites) for sites in filtered_sites.values())}")
    print(f"After specificity filter: {sum(len(sites) for sites in final_sites.values())}")


if __name__ == '__main__':
    main()
