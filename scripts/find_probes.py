import os
import sys
import argparse
import shutil
import logging
from datetime import datetime
from typing import List, Dict
import time

# Add src directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from probe_designer.config import ConfigManager
from probe_designer.database import DatabaseInterface
from probe_designer.search_strategies import BindingSiteSearcher
from probe_designer.filtering import SequenceFilter
from probe_designer.utils import load_gene_list


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
    parser = argparse.ArgumentParser(description="Find binding sites using brute force strategy (works with Ensembl and NCBI)")
    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument('--genes_file', help='Path to a text file of gene names (one per line)')
    gene_group.add_argument('--gene_info', help='Excel file with gene_name, GI columns (auto-detects NCBI IDs)')
    parser.add_argument('--config', default='configs/config_bruteforce.yaml', help='Path to configuration file')
    parser.add_argument('--species', default='mouse', help='Species name (mouse, human, rat, zebrafish)')
    parser.add_argument('--database', choices=['ensembl', 'ncbi'], help='Database type (overrides config)')
    parser.add_argument('--output_dir', help='Override output directory from config')
    parser.add_argument('--genome_fasta', help='Override genome FASTA path from config')
    parser.add_argument('--blast_species', nargs='+', help='Target species for BLAST search (default: from config blast.species)')
    args = parser.parse_args()

    # Load configuration
    cfg = ConfigManager(args.config, species=args.species)
    
    # Override config with command line arguments if provided
    if args.database:
        cfg.database.database_type = args.database
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

    # Prepare genome accessor using priority chain
    accessor = None
    try:
        accessor = db.initialize_genome_accessor(cfg.genome)
    except RuntimeError as e:
        print(f"Warning: genome accessor not available ({e}), proceeding without it")

    # Read gene list (from --genes_file or --gene_info)
    gene_ids = None
    if args.gene_info:
        import pandas as pd
        df = pd.read_excel(args.gene_info)
        genes = df['gene_name'].astype(str).str.strip().tolist()
        if 'GI' in df.columns:
            gene_ids = df['GI'].astype(str).str.strip().tolist()
            cfg.database.database_type = 'ncbi'
    else:
        genes = load_gene_list(args.genes_file)
    print(f"Processing {len(genes)} genes using {cfg.database.database_type} database...")

    # Get gene sequences (unpack result structure)
    seq_result = db.get_gene_sequences(genes, gene_ids=gene_ids)
    sequences = seq_result.get('sequences', {}) if isinstance(seq_result, dict) else seq_result
    errors = seq_result.get('errors', {}) if isinstance(seq_result, dict) else {}
    print(f"Retrieved sequences for {len(sequences)} genes")
    # Debug: show brief sequence info
    try:
        example_keys = list(sequences.keys())[:3]
        print(f"Example genes: {example_keys}")
        for g in example_keys:
            seq_obj = sequences.get(g, {}) or {}
            seq_str = seq_obj.get('sequence', '')
            print(f"  {g}: length={len(seq_str)} preview={seq_str[:30]}...")
    except Exception as e:
        print(f"Debug print of sequences failed: {e}")
    # Quick diagnostics for missing sequences
    missing = [g for g in genes if g not in sequences or not sequences[g].get('sequence')]
    if missing:
        print(f"Missing sequences for {len(missing)} genes: {missing[:5]}{'...' if len(missing) > 5 else ''}")
    # Show first error per gene if any
    if errors:
        sample_errors = {g: (errs[0] if errs else '') for g, errs in errors.items() if errs}
        if sample_errors:
            print("Sample errors:")
            for g, msg in list(sample_errors.items())[:5]:
                print(f"  {g}: {msg}")

    # Initialize searcher
    searcher = BindingSiteSearcher(cfg.search, cfg.filter, genome_accessor=accessor)
    print("Starting binding site search (bruteforce)...")
    t0 = time.time()
    binding_sites = searcher.search_all_genes(sequences)
    print(f"Binding site search completed in {time.time() - t0:.1f}s")

    # Create output directory
    os.makedirs(cfg.output.output_dir, exist_ok=True)

    # Enrich sites with padlock arms and target sequence expected by filters
    enriched_sites: Dict[str, List[Dict]] = {}
    bds_len = cfg.search.binding_site_length
    for gene, sites in binding_sites.items():
        enriched_sites[gene] = []
        target_seq = sequences.get(gene, {}).get('sequence', '')
        for site in sites:
            probe_seq = site.get('sequence', '')
            arm_3prime = probe_seq[: bds_len // 2]  # front half
            arm_5prime = probe_seq[bds_len // 2 :]  # back half
            new_site = dict(site)
            new_site['arm_3prime'] = arm_3prime
            new_site['arm_5prime'] = arm_5prime
            new_site['target_sequence'] = probe_seq
            enriched_sites[gene].append(new_site)

    # Initialize filter and run full pipeline (batch APIs)
    flt = SequenceFilter(cfg.filter, cfg.blast)
    print("\nStarting BLAST analysis...")

    # 1) Pre-BLAST thermal filtering (batch)
    print("Applying pre-BLAST filters...")
    filtered_sites = flt.pre_blast_filter(enriched_sites)

    # 2) Run BLAST
    print("Running BLAST analysis...")
    blast_results_path = flt.run_blast(filtered_sites, cfg.output.output_dir)

    # 3) Specificity filtering (batch)
    print("Applying specificity filtering...")
    final_sites = flt.specificity_filter(filtered_sites, blast_results_path)

    # 4) Save results
    print("Saving final results...")
    flt.save_filtered_results(final_sites, cfg.output.output_dir)
    # Also save original sites for comparison
    searcher.save_binding_sites(binding_sites, cfg.output.output_dir)

    print(f"Results saved to {cfg.output.output_dir}")
    print(f"Original sites: {sum(len(sites) for sites in binding_sites.values())}")
    print(f"After pre-BLAST filter: {sum(len(sites) for sites in filtered_sites.values())}")
    print(f"After specificity filter: {sum(len(sites) for sites in final_sites.values())}")


if __name__ == '__main__':
    main()
