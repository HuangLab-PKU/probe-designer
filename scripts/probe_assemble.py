#!/usr/bin/env python3
"""
Probe Assembly Script

This script assembles probes from binding sites and backbone sequences.
Uses the simplified ProbeAssembler class.

Author: Refactored version
"""

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List, Any

# Add code directory to path
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from probe_designer.config import ProbeConfig
from probe_designer.probe_assembly import ProbeAssembler


def load_binding_sites_from_json(json_file: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Load binding sites data from JSON file
    
    Args:
        json_file: Path to JSON file
        
    Returns:
        Dict mapping gene_name to list of binding site dicts
    """
    if not os.path.exists(json_file):
        raise FileNotFoundError(f"Binding sites file not found: {json_file}")
    
    with open(json_file, 'r', encoding='utf-8') as f:
        binding_sites = json.load(f)
    
    print(f"Loaded {sum(len(sites) for sites in binding_sites.values())} binding sites from {len(binding_sites)} genes")
    return binding_sites


def load_binding_sites_from_excel(excel_file: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Load binding sites data from Excel file (filtered results)
    
    Args:
        excel_file: Path to Excel file
        
    Returns:
        Dict mapping gene_name to list of binding site dicts
    """
    import pandas as pd
    
    if not os.path.exists(excel_file):
        raise FileNotFoundError(f"Binding sites file not found: {excel_file}")
    
    df = pd.read_excel(excel_file)
    
    # Validate required columns
    required_cols = ['gene_name', 'arm3', 'arm5', 'st', 'en', 'g_content', 'tm', 'tm3', 'tm5']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Excel file missing required columns: {missing_cols}")
    
    # Convert to required format
    binding_sites = {}
    for _, row in df.iterrows():
        gene_name = str(row['gene_name'])
        if gene_name not in binding_sites:
            binding_sites[gene_name] = []
        
        site = {
            'arm_3prime': str(row['arm3']),
            'arm_5prime': str(row['arm5']),
            'st': int(row['st']),
            'en': int(row['en']),
            'g_content': float(row['g_content']),
            'tm': float(row['tm']),
            'tm_3prime': float(row['tm3']),
            'tm_5prime': float(row['tm5'])
        }
        
        # Optional fields
        if 'isoform_overlap_num' in row:
            site['isoform_overlap_num'] = row['isoform_overlap_num']
            
        binding_sites[gene_name].append(site)
    
    print(f"Loaded {sum(len(sites) for sites in binding_sites.values())} binding sites from {len(binding_sites)} genes")
    return binding_sites


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Probe assembly tool")
    parser.add_argument("--binding_sites", required=True, help="Binding sites file path (JSON or Excel)")
    parser.add_argument("--gene_info", required=True, help="Gene info Excel file path (must contain gene_name and No. columns)")
    parser.add_argument("--backbone_file", required=True, help="Backbone Excel file path (must contain No. and Sequence columns)")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--ilock", choices=["3prime", "5prime", "both"], help="Add iLock modification (optional)")
    
    args = parser.parse_args()
    
    # Determine output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Output directory: {output_dir}")
    
    # Load binding sites
    binding_sites_file = args.binding_sites
    if binding_sites_file.endswith('.json'):
        binding_sites = load_binding_sites_from_json(binding_sites_file)
    elif binding_sites_file.endswith(('.xlsx', '.xls')):
        binding_sites = load_binding_sites_from_excel(binding_sites_file)
    else:
        raise ValueError(f"Unsupported file format: {binding_sites_file}. Please use JSON or Excel files.")
    
    if not binding_sites:
        print("Warning: No binding sites found")
        return
    
    # Create ProbeConfig
    probe_config = ProbeConfig(backbone_file=args.backbone_file)
    
    # Create ProbeAssembler
    assembler = ProbeAssembler(probe_config)
    
    # Assemble probes
    print("Starting probe assembly...")
    probe_df = assembler.assemble_probes(binding_sites, args.gene_info)
    probe_df['backbone_file'] = os.path.basename(args.backbone_file)
    if probe_df.empty:
        print("Warning: No probes generated")
        return
    # Optional: Add iLock modification
    if args.ilock:
        print(f"Adding iLock modification (position: {args.ilock})...")
        probe_df = assembler.add_ilock_modification(probe_df, position=args.ilock)
    
    # Save probes
    print("Saving probe results...")
    assembler.save_probes(probe_df, output_dir)
    
    print(f"Probe assembly completed! Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
