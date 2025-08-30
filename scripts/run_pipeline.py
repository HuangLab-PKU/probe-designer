#!/usr/bin/env python3
"""
Multi-database probe design pipeline
Runs probe design with multiple databases and merges results

Usage:
    python run_pipeline.py --genes-file genes.txt --project-id project_001 --strategy consensus
    python run_pipeline.py --genes-file genes.txt --project-id project_002 --strategy specific --databases ensembl ncbi
"""

import os
import sys
import argparse
import subprocess
import json
import pandas as pd
import shutil
from pathlib import Path
from typing import List, Dict, Optional

# Add src directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.config import ConfigManager


def read_gene_list(path: str) -> List[str]:
    """Read gene list from file."""
    with open(path, 'r', encoding='utf-8') as f:
        return [line.strip() for line in f if line.strip()]


def run_probe_design(genes_file: str, config_file: str, output_dir: str, 
                    database: str, strategy: str, blast_species: Optional[List[str]] = None) -> bool:
    """Run probe design with specified parameters."""
    
    # Determine which script to run based on strategy
    if strategy == "consensus":
        script = "find_consensus_probes.py"
        default_config = "configs/config_consensus.yaml"
    elif strategy == "specific":
        script = "find_specific_probes.py"
        default_config = "configs/config_specific.yaml"
    elif strategy == "bruteforce":
        script = "find_probes.py"
        default_config = "configs/config_bruteforce.yaml"
    else:
        print(f"Unknown strategy: {strategy}")
        return False
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Resolve which config will be used
    used_config = config_file or default_config

    # Build command
    cmd = [
        sys.executable, f"scripts/{script}",
        "--genes_file", genes_file,
        "--config", used_config,
        "--output_dir", output_dir
    ]
    
    # Add database override for bruteforce strategy
    if strategy == "bruteforce":
        cmd.extend(["--database", database])
    
    # Add blast species if provided
    if blast_species:
        cmd.extend(["--blast_species"] + blast_species)
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        # Stream output in real-time for better progress visibility
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
        process.wait()
        if process.returncode != 0:
            print(f"Command failed with exit code {process.returncode}")
            return False

        # Copy the used config into the output directory for record keeping
        try:
            if os.path.exists(used_config):
                dst_config = os.path.join(output_dir, "config_used.yaml")
                shutil.copy2(used_config, dst_config)
                print(f"Saved used config to: {dst_config}")
        except Exception as copy_err:
            print(f"Warning: failed to copy config file: {copy_err}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print("Output:")
        print(e.stdout)
        print("Errors:")
        print(e.stderr)
        return False


def check_missing_genes(merged_results_file: str, genes_file: str, organism: str = "mouse") -> List[str]:
    """Check for missing genes using merge_results.py."""
    
    cmd = [
        sys.executable, "test/merge_results.py",
        "--results-dir", os.path.dirname(merged_results_file),
        "--gene-list", genes_file,
        "--organism", organism,
        "--output", merged_results_file,
        "--missing-output", "missing_genes.txt"
    ]
    
    print(f"Checking missing genes: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Missing genes check output:")
        print(result.stdout)
        
        # Read missing genes file if it exists
        missing_genes = []
        if os.path.exists("missing_genes.txt"):
            with open("missing_genes.txt", 'r') as f:
                missing_genes = [line.strip() for line in f if line.strip()]
        
        return missing_genes
        
    except subprocess.CalledProcessError as e:
        print(f"Missing genes check failed: {e}")
        return []


def merge_results(results_dirs: List[str], output_file: str) -> bool:
    """Merge results from multiple directories."""
    
    # Create a temporary directory for merging
    temp_dir = "temp_merge"
    os.makedirs(temp_dir, exist_ok=True)
    
    # Copy all results to temp directory
    for i, results_dir in enumerate(results_dirs):
        if os.path.exists(results_dir):
            # Copy files to temp directory with prefix
            for file in os.listdir(results_dir):
                if file.endswith(('.json', '.fasta', '.xlsx')):
                    src = os.path.join(results_dir, file)
                    dst = os.path.join(temp_dir, f"db_{i}_{file}")
                    if os.path.isfile(src):
                        import shutil
                        shutil.copy2(src, dst)
    
    # Run merge_results.py
    cmd = [
        sys.executable, "test/merge_results.py",
        "--results-dir", temp_dir,
        "--output", output_file
    ]
    
    print(f"Merging results: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Merge output:")
        print(result.stdout)
        
        # Clean up temp directory
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Merge failed: {e}")
        # Clean up temp directory
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Multi-database probe design pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run consensus strategy with Ensembl only
  python run_pipeline.py --genes-file genes.txt --project-id project_001 --strategy consensus
  
  # Run specific strategy with both databases
  python run_pipeline.py --genes-file genes.txt --project-id project_002 --strategy specific --databases ensembl ncbi
  
  # Run bruteforce strategy with custom config
  python run_pipeline.py --genes-file genes.txt --project-id project_003 --strategy bruteforce --databases ensembl ncbi --config my_config.yaml
        """
    )
    
    parser.add_argument("--genes-file", required=True, help="Path to gene list file")
    parser.add_argument("--project-id", required=True, help="Project identifier")
    parser.add_argument("--strategy", choices=["consensus", "specific", "bruteforce"], 
                       default="consensus", help="Probe design strategy")
    parser.add_argument("--databases", nargs="+", choices=["ensembl", "ncbi"], 
                       default=["ensembl"], help="Databases to use")
    parser.add_argument("--config", help="Custom configuration file")
    parser.add_argument("--blast-species", nargs="+", help="Target species for BLAST search")
    parser.add_argument("--organism", default="mouse", choices=["mouse", "human", "rat", "zebrafish"],
                       help="Target organism")
    parser.add_argument("--skip-merge", action="store_true", help="Skip result merging")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.genes_file):
        print(f"Error: Gene list file not found: {args.genes_file}")
        sys.exit(1)
    
    # Create project directory
    project_dir = f"results/{args.project_id}"
    os.makedirs(project_dir, exist_ok=True)
    
    print(f"Starting pipeline for project: {args.project_id}")
    print(f"Strategy: {args.strategy}")
    print(f"Databases: {args.databases}")
    print(f"Organism: {args.organism}")
    
    # Run probe design for each database
    results_dirs = []
    for i, database in enumerate(args.databases):
        print(f"\n{'='*50}")
        print(f"Running {database} database ({i+1}/{len(args.databases)})")
        print(f"{'='*50}")
        
        # Create database-specific output directory
        db_output_dir = os.path.join(project_dir, f"{database}_results")
        
        # Run probe design
        success = run_probe_design(
            genes_file=args.genes_file,
            config_file=args.config,
            output_dir=db_output_dir,
            database=database,
            strategy=args.strategy,
            blast_species=args.blast_species
        )
        
        if success:
            results_dirs.append(db_output_dir)
            print(f"✓ {database} database completed successfully")
        else:
            print(f"✗ {database} database failed")
    
    if not results_dirs:
        print("No successful runs. Exiting.")
        sys.exit(1)
    
    # Merge results if multiple databases were used
    if len(results_dirs) > 1 and not args.skip_merge:
        print(f"\n{'='*50}")
        print("Merging results from multiple databases")
        print(f"{'='*50}")
        
        merged_file = os.path.join(project_dir, "merged_results.xlsx")
        merge_success = merge_results(results_dirs, merged_file)
        
        if merge_success:
            print(f"✓ Results merged successfully: {merged_file}")
            
            # Check for missing genes
            print(f"\n{'='*50}")
            print("Checking for missing genes")
            print(f"{'='*50}")
            
            missing_genes = check_missing_genes(merged_file, args.genes_file, args.organism)
            
            if missing_genes:
                print(f"Found {len(missing_genes)} missing genes:")
                for gene in missing_genes:
                    print(f"  - {gene}")
                
                # Save missing genes for potential re-run
                missing_file = os.path.join(project_dir, "missing_genes.txt")
                with open(missing_file, 'w') as f:
                    for gene in missing_genes:
                        f.write(f"{gene}\n")
                print(f"Missing genes saved to: {missing_file}")
            else:
                print("✓ All genes have probes designed!")
        else:
            print("✗ Result merging failed")
    
    # Generate summary
    print(f"\n{'='*50}")
    print("Pipeline Summary")
    print(f"{'='*50}")
    print(f"Project ID: {args.project_id}")
    print(f"Strategy: {args.strategy}")
    print(f"Databases run: {args.databases}")
    print(f"Successful runs: {len(results_dirs)}")
    print(f"Results directory: {project_dir}")
    
    if len(results_dirs) > 1 and not args.skip_merge:
        print(f"Merged results: {os.path.join(project_dir, 'merged_results.xlsx')}")
    
    print(f"\nPipeline completed successfully!")


if __name__ == "__main__":
    main()
