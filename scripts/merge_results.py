#!/usr/bin/env python3
"""
Merge Results Script

Traverse subdirectories under results directory, collect all filtered_binding_sites.xlsx files,
filter, deduplicate, and sort according to gene_info, then output merged results.

Usage:
    python merge_results.py --results-dir results/ --gene-info gene_info.xlsx --output merged_binding_sites.xlsx
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from typing import List, Set, Dict

# Add src directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def find_filtered_binding_sites_files(results_dir: str) -> List[str]:
    """
    Traverse all subdirectories under results directory and find filtered_binding_sites.xlsx files
    
    Args:
        results_dir: Results directory path
        
    Returns:
        List of all found filtered_binding_sites.xlsx file paths
    """
    binding_sites_files = []
    results_path = Path(results_dir)
    
    if not results_path.exists():
        raise FileNotFoundError(f"Results directory does not exist: {results_dir}")
    
    # Traverse all subdirectories
    for subdir in results_path.iterdir():
        if subdir.is_dir():
            # Look for filtered_binding_sites.xlsx in each subdirectory
            binding_sites_file = subdir / "filtered_binding_sites.xlsx"
            if binding_sites_file.exists():
                binding_sites_files.append(str(binding_sites_file))
                print(f"Found: {binding_sites_file}")
    
    return binding_sites_files


def load_gene_info(gene_info_file: str) -> pd.DataFrame:
    """
    Load gene_info Excel file
    
    Args:
        gene_info_file: Path to gene_info Excel file (must contain gene_name and No. columns)
        
    Returns:
        DataFrame with gene_name and No. columns
    """
    if not os.path.exists(gene_info_file):
        raise FileNotFoundError(f"Gene info file not found: {gene_info_file}")
    
    df = pd.read_excel(gene_info_file)
    
    # Validate required columns
    if 'gene_name' not in df.columns or 'No.' not in df.columns:
        raise ValueError(f"Gene info file must contain 'gene_name' and 'No.' columns. Found: {df.columns.tolist()}")
    
    # Normalize gene names: strip whitespace
    df['gene_name'] = df['gene_name'].astype(str).str.strip()
    
    return df


def load_and_merge_binding_sites(binding_sites_files: List[str]) -> pd.DataFrame:
    """
    Load and merge all binding sites files
    
    Args:
        binding_sites_files: List of binding sites Excel file paths
        
    Returns:
        Merged DataFrame
    """
    if not binding_sites_files:
        raise ValueError("No binding sites files found")
    
    all_dfs = []
    
    for file_path in binding_sites_files:
        try:
            df = pd.read_excel(file_path)
            # Validate required columns
            required_cols = ['gene_name']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                print(f"Warning: {file_path} missing required columns: {missing_cols}, skipping this file")
                continue
            
            # Normalize gene names
            df['gene_name'] = df['gene_name'].astype(str).str.strip()
            all_dfs.append(df)
            print(f"Loaded: {file_path} ({len(df)} rows)")
        except Exception as e:
            print(f"Warning: Error loading {file_path}: {e}, skipping this file")
            continue
    
    if not all_dfs:
        raise ValueError("No valid binding sites files could be loaded")
    
    # Merge all DataFrames
    merged_df = pd.concat(all_dfs, ignore_index=True)
    print(f"Merged total: {len(merged_df)} rows from {len(binding_sites_files)} files")
    
    return merged_df


def filter_by_gene_info(df: pd.DataFrame, gene_info: pd.DataFrame) -> pd.DataFrame:
    """
    Filter binding sites by gene_info, keeping only genes present in gene_info
    
    Args:
        df: Merged binding sites DataFrame
        gene_info: gene_info DataFrame
        
    Returns:
        Filtered DataFrame
    """
    # Get gene list from gene_info (case-insensitive matching)
    gene_info_genes = set(gene_info['gene_name'].str.upper())
    
    # Create case-insensitive matching map
    gene_name_upper = df['gene_name'].str.upper()
    mask = gene_name_upper.isin(gene_info_genes)
    
    filtered_df = df[mask].copy()
    
    print(f"Filtered by gene_info: {len(df)} -> {len(filtered_df)} rows")
    
    return filtered_df


def deduplicate_binding_sites(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate binding sites (based on gene_name, st, en)
    
    Args:
        df: binding sites DataFrame
        
    Returns:
        Deduplicated DataFrame
    """
    if 'st' in df.columns and 'en' in df.columns:
        # Deduplicate based on gene_name, st, en
        before_count = len(df)
        df_dedup = df.drop_duplicates(subset=['gene_name', 'st', 'en'], keep='first')
        after_count = len(df_dedup)
        print(f"Deduplication: {before_count} -> {after_count} rows (based on gene_name, st, en)")
    else:
        # If no st, en, deduplicate based on gene_name, sequence
        if 'sequence' in df.columns:
            before_count = len(df)
            df_dedup = df.drop_duplicates(subset=['gene_name', 'sequence'], keep='first')
            after_count = len(df_dedup)
            print(f"Deduplication: {before_count} -> {after_count} rows (based on gene_name, sequence)")
        else:
            print("Warning: Cannot deduplicate, missing st/en or sequence columns")
            df_dedup = df
    
    return df_dedup


def sort_by_gene_info(df: pd.DataFrame, gene_info: pd.DataFrame) -> pd.DataFrame:
    """
    Sort by gene_info order (by row order in gene_info, not by No.)
    Then sort by st and en
    
    Args:
        df: binding sites DataFrame
        gene_info: gene_info DataFrame
        
    Returns:
        Sorted DataFrame
    """
    # Create copy to avoid SettingWithCopyWarning
    df = df.copy()
    
    # Create gene_name -> index position in gene_info mapping (case-insensitive)
    gene_to_order = {}
    for idx, row in gene_info.iterrows():
        gene_name = str(row['gene_name']).strip()
        gene_to_order[gene_name.upper()] = idx  # Use original index as sort key
    
    # Add sort key column (by order in gene_info)
    df['_sort_key'] = df['gene_name'].str.upper().map(gene_to_order)
    
    # For genes not matched, use a large number to put them at the end
    df['_sort_key'] = df['_sort_key'].fillna(999999)
    
    # Sort: first by gene_info order, then by st, finally by en
    sort_columns = ['_sort_key']
    if 'st' in df.columns:
        sort_columns.append('st')
    if 'en' in df.columns:
        sort_columns.append('en')
    
    df_sorted = df.sort_values(by=sort_columns).drop(columns=['_sort_key'])
    
    print(f"Sorted by gene_info order (then by st, en)")
    
    return df_sorted


def find_missing_genes(gene_info: pd.DataFrame, merged_df: pd.DataFrame) -> List[str]:
    """
    Find genes present in gene_info but not found in merged results
    Return in gene_info order (no re-sorting)
    
    Args:
        gene_info: gene_info DataFrame
        merged_df: Merged binding sites DataFrame
        
    Returns:
        List of missing genes (in gene_info order)
    """
    merged_genes = set(merged_df['gene_name'].str.upper())
    
    # Traverse in gene_info order to find missing genes
    missing_genes = []
    for _, row in gene_info.iterrows():
        gene_name = str(row['gene_name']).strip()
        if gene_name.upper() not in merged_genes:
            missing_genes.append(gene_name)  # Use original case
    
    return missing_genes


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Merge binding sites files from multiple result directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python merge_results.py --results-dir results/ --gene-info gene_info.xlsx --output merged_binding_sites.xlsx
  python merge_results.py --results-dir results/ --gene-info gene_info.xlsx --output merged.xlsx --missing-output missing_genes.txt
        """
    )
    
    parser.add_argument("--results-dir", "-r", required=True, 
                       help="Results directory path (will traverse all subdirectories)")
    parser.add_argument("--gene-info", "-g", required=True,
                       help="Gene info Excel file path (must contain gene_name and No. columns)")
    parser.add_argument("--output", "-o", required=True,
                       help="Output file path (merged binding sites Excel file)")
    parser.add_argument("--missing-output", "-m",
                       help="Missing genes output file path (default: same directory as output file)")
    
    args = parser.parse_args()
    
    try:
        # 1. Find all filtered_binding_sites.xlsx files
        print(f"Searching for filtered_binding_sites.xlsx files in {args.results_dir}...")
        binding_sites_files = find_filtered_binding_sites_files(args.results_dir)
        
        if not binding_sites_files:
            print(f"Error: No filtered_binding_sites.xlsx files found in {args.results_dir}")
            sys.exit(1)
        
        print(f"Found {len(binding_sites_files)} files\n")
        
        # 2. Load and merge all binding sites
        print("Loading and merging binding sites...")
        merged_df = load_and_merge_binding_sites(binding_sites_files)
        
        # 3. Load gene_info
        print(f"\nLoading gene_info: {args.gene_info}")
        gene_info = load_gene_info(args.gene_info)
        print(f"Gene info contains {len(gene_info)} genes\n")
        
        # 4. Filter by gene_info
        print("Filtering by gene_info...")
        filtered_df = filter_by_gene_info(merged_df, gene_info)
        
        # 5. Deduplicate
        print("\nDeduplicating...")
        dedup_df = deduplicate_binding_sites(filtered_df)
        
        # 6. Sort by gene_info
        print("\nSorting...")
        sorted_df = sort_by_gene_info(dedup_df, gene_info)
        
        # 7. Determine output file path
        output_path = Path(args.output)
        
        # If output path is a directory or has no .xlsx/.xls extension, auto-add filename
        if output_path.is_dir() or (not output_path.suffix) or output_path.suffix.lower() not in ['.xlsx', '.xls']:
            if output_path.is_dir() or (not output_path.suffix):
                # If directory or no extension, add filename
                output_path = output_path / "merged_binding_sites.xlsx"
            else:
                # If has extension but not Excel format, replace with .xlsx
                output_path = output_path.with_suffix('.xlsx')
        
        # Ensure parent directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save results
        sorted_df.to_excel(str(output_path), index=False)
        print(f"\nMerged results saved to: {output_path}")
        print(f"  Total: {len(sorted_df)} rows, {sorted_df['gene_name'].nunique()} genes")
        
        # 8. Find missing genes
        print("\nChecking for missing genes...")
        missing_genes = find_missing_genes(gene_info, sorted_df)
        
        if missing_genes:
            # Determine missing genes output file path
            if args.missing_output:
                missing_output_path = args.missing_output
            else:
                # Default: same directory as output file
                missing_output_path = str(output_path.parent / "missing_genes.txt")
            
            # Save missing genes list
            with open(missing_output_path, 'w', encoding='utf-8') as f:
                for gene in missing_genes:
                    f.write(f"{gene}\n")
            
            print(f"Found {len(missing_genes)} missing genes:")
            for gene in missing_genes[:20]:  # Show first 20
                print(f"  - {gene}")
            if len(missing_genes) > 20:
                print(f"  ... (total {len(missing_genes)} genes)")
            print(f"Missing genes list saved to: {missing_output_path}")
        else:
            print("✓ All genes have binding sites!")
        
        # 9. Display statistics
        print(f"\nStatistics:")
        print(f"  Total binding sites: {len(sorted_df)}")
        print(f"  Number of genes: {sorted_df['gene_name'].nunique()}")
        if 'g_content' in sorted_df.columns:
            print(f"  Average G content: {sorted_df['g_content'].mean():.2f}")
        if 'tm' in sorted_df.columns:
            print(f"  Average Tm: {sorted_df['tm'].mean():.1f} °C")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
