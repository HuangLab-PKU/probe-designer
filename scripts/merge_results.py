#!/usr/bin/env python3
"""
Results merging script
Merge probe design results from multiple result directories

Usage:
    python merge_results.py --results-dir results/ --gene-list gene_list.xlsx --output merged_results.xlsx
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path

# Add src directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from utils import (
    load_gene_list, merge_results_from_directories, 
    find_missing_genes, save_missing_genes, format_gene_name
)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Merge probe design results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python merge_results.py --results-dir results/ --gene-list genes.xlsx --output merged.xlsx
  python merge_results.py --results-dir results/ --organism mouse --output merged.xlsx
        """
    )
    
    parser.add_argument("--results-dir", "-r", required=True,
                       help="Results directory path")
    parser.add_argument("--gene-list", "-g",
                       help="Original gene list file")
    parser.add_argument("--organism", "-o", 
                       choices=["mouse", "human"], default="mouse",
                       help="Species name (default: mouse)")
    parser.add_argument("--output", "-out", required=True,
                       help="Output file path")
    parser.add_argument("--missing-output",
                       help="Missing genes output file path")
    
    args = parser.parse_args()
    
    try:
        # Check results directory
        if not os.path.exists(args.results_dir):
            print(f"Error: Results directory does not exist: {args.results_dir}")
            sys.exit(1)
        
        # Merge results
        print(f"Merging results directory: {args.results_dir}")
        merged_df = merge_results_from_directories(args.results_dir, args.output)
        
        print(f"Successfully merged {len(merged_df)} probes")
        print(f"Results saved to: {args.output}")
        
        # If gene list is provided, check for missing genes
        if args.gene_list:
            print(f"Checking missing genes: {args.gene_list}")
            
            # Load original gene list
            gene_list = load_gene_list(args.gene_list)
            
            # Format gene names
            formatted_genes = []
            for gene in gene_list:
                formatted_gene = format_gene_name(gene, args.organism)
                formatted_genes.append(formatted_gene)
            
            # Find missing genes
            missing_genes = find_missing_genes(formatted_genes, merged_df)
            
            if missing_genes:
                print(f"Found {len(missing_genes)} missing genes:")
                for gene in missing_genes:
                    print(f"  - {gene}")
                
                # Save missing genes list
                missing_output = args.missing_output or "missing_genes.txt"
                save_missing_genes(missing_genes, missing_output)
                print(f"Missing genes list saved to: {missing_output}")
            else:
                print("All genes successfully designed probes!")
        
        # Display statistics
        print(f"\nStatistics:")
        print(f"  Total probes: {len(merged_df)}")
        print(f"  Genes involved: {merged_df['gene_name'].nunique()}")
        print(f"  Average probes per gene: {len(merged_df) / merged_df['gene_name'].nunique():.1f}")
        
        if 'g_content' in merged_df.columns:
            print(f"  Average G content: {merged_df['g_content'].mean():.2f}")
        if 'tm' in merged_df.columns:
            print(f"  Average melting temperature: {merged_df['tm'].mean():.1f} °C")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
