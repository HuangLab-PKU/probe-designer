#!/usr/bin/env python3
"""
Main entry point for the automated DNA probe design software.
"""

import argparse
import sys
import os
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from config import ConfigManager
from database import DatabaseInterface
from search_strategies import BindingSiteSearcher
from filtering import SequenceFilter
from probe_assembly import ProbeAssembler, ProbeValidator
from utils import (
    create_output_directory, load_gene_list, save_gene_list,
    log_message, check_dependencies, print_dependency_status,
    create_config_template, calculate_statistics, save_statistics
)


class ProbeDesigner:
    """Main orchestrator class for the probe design workflow."""
    
    def __init__(self, config_file: str = None):
        self.config_manager = ConfigManager(config_file)
        self.config = self.config_manager.config
        
        # Initialize components
        self.database = DatabaseInterface(self.config.database)
        self.searcher = BindingSiteSearcher(self.config.search)
        self.filter = SequenceFilter(self.config.filter, self.config.blast)
        self.assembler = ProbeAssembler(self.config.probe)
        self.validator = ProbeValidator()
        
        # Setup output directory
        self.output_dir = create_output_directory(
            self.config.output.output_dir,
            self.config.output.create_timestamp
        )
        self.log_file = os.path.join(self.output_dir, "probe_design.log")
    
    def run(self, gene_list_file: str):
        """Execute the complete probe design workflow."""
        log_message("Starting probe design workflow", self.log_file)
        
        # Load gene list
        log_message(f"Loading gene list from: {gene_list_file}", self.log_file)
        gene_list = load_gene_list(gene_list_file, "gene_name")
        log_message(f"Loaded {len(gene_list)} genes", self.log_file)
        
        # Save processed gene list
        processed_gene_file = os.path.join(self.output_dir, "gene_name_list_tosearch.txt")
        save_gene_list(gene_list, processed_gene_file)
        
        # Fetch sequences
        log_message("Fetching gene sequences from database", self.log_file)
        sequences = self.database.fetch_sequences(gene_list)
        log_message(f"Retrieved sequences for {len(sequences)} genes", self.log_file)
        
        # Search for binding sites
        log_message("Searching for binding sites", self.log_file)
        binding_sites = self.searcher.search_all_genes(sequences)
        log_message(f"Found binding sites for {len(binding_sites)} genes", self.log_file)
        
        # Pre-BLAST filtering
        log_message("Applying pre-BLAST filters", self.log_file)
        filtered_sites = self.filter.pre_blast_filter(binding_sites)
        log_message(f"Pre-filtered to {sum(len(sites) for sites in filtered_sites.values())} sites", self.log_file)
        
        # Run BLAST
        log_message("Running BLAST analysis", self.log_file)
        blast_results_file = self.filter.run_blast(filtered_sites, self.output_dir)
        log_message(f"BLAST completed: {blast_results_file}", self.log_file)
        
        # Post-BLAST filtering
        log_message("Applying post-BLAST filters", self.log_file)
        final_sites = self.filter.post_blast_filter(filtered_sites, blast_results_file)
        log_message(f"Post-filtered to {sum(len(sites) for sites in final_sites.values())} sites", self.log_file)
        
        # Save filtered results
        log_message("Saving filtered results", self.log_file)
        self.filter.save_filtered_results(final_sites, self.output_dir)
        
        # Assemble probes
        log_message("Assembling probes", self.log_file)
        probe_df = self.assembler.assemble_probes(final_sites)
        log_message(f"Assembled {len(probe_df)} probes", self.log_file)
        
        # Save probes
        log_message("Saving probe results", self.log_file)
        self.assembler.save_probes(probe_df, self.output_dir)
        
        # Validate probes
        log_message("Validating probes", self.log_file)
        validation_result = self.validator.validate_probe_set(probe_df)
        if not validation_result['is_valid']:
            log_message("WARNING: Probe validation failed", self.log_file)
            for error in validation_result['errors']:
                log_message(f"  Error: {error}", self.log_file)
        
        # Calculate and save statistics
        log_message("Calculating statistics", self.log_file)
        stats = calculate_statistics(final_sites)
        stats_file = os.path.join(self.output_dir, "statistics.json")
        save_statistics(stats, stats_file)
        
        log_message("Probe design workflow completed successfully", self.log_file)
        log_message(f"Results saved to: {self.output_dir}", self.log_file)
        
        return probe_df


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Automated DNA probe design software",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s gene_list.xlsx
  %(prog)s genes.txt --config my_config.json
  %(prog)s --create-config template.json
  %(prog)s --check-deps
        """
    )
    
    parser.add_argument(
        "gene_list_file",
        nargs="?",
        help="Path to gene list file (.xlsx, .csv, or .txt)"
    )
    
    parser.add_argument(
        "--config", "-c",
        help="Path to configuration JSON file"
    )
    
    parser.add_argument(
        "--create-config",
        help="Create a configuration template file"
    )
    
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help="Check if all required dependencies are installed"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 1.0.0"
    )
    
    args = parser.parse_args()
    
    # Handle special commands
    if args.check_deps:
        print_dependency_status()
        return
    
    if args.create_config:
        create_config_template(args.create_config)
        return
    
    # Validate required arguments
    if not args.gene_list_file:
        parser.error("gene_list_file is required")
    
    if not os.path.exists(args.gene_list_file):
        parser.error(f"Gene list file not found: {args.gene_list_file}")
    
    # Check dependencies
    deps = check_dependencies()
    missing = [pkg for pkg, installed in deps.items() if not installed]
    if missing:
        print("ERROR: Missing required dependencies:")
        for pkg in missing:
            print(f"  - {pkg}")
        print("\nInstall missing packages with:")
        print(f"pip install {' '.join(missing)}")
        sys.exit(1)
    
    try:
        # Run probe design
        designer = ProbeDesigner(args.config)
        probe_df = designer.run(args.gene_list_file)
        
        print(f"\nProbe design completed successfully!")
        print(f"Results saved to: {designer.output_dir}")
        print(f"Total probes generated: {len(probe_df)}")
        
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
