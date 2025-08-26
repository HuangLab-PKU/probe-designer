#!/usr/bin/env python3
"""
结果合并脚本
合并多个结果目录中的探针设计结果

使用方法:
    python merge_results.py --results-dir results/ --gene-list gene_list.xlsx --output merged_results.xlsx
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path

# 添加src目录到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from utils import (
    load_gene_list, merge_results_from_directories, 
    find_missing_genes, save_missing_genes, format_gene_name
)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="合并探针设计结果",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python merge_results.py --results-dir results/ --gene-list genes.xlsx --output merged.xlsx
  python merge_results.py --results-dir results/ --organism mouse --output merged.xlsx
        """
    )
    
    parser.add_argument("--results-dir", "-r", required=True,
                       help="结果目录路径")
    parser.add_argument("--gene-list", "-g",
                       help="原始基因列表文件")
    parser.add_argument("--organism", "-o", 
                       choices=["mouse", "human"], default="mouse",
                       help="物种名称 (默认: mouse)")
    parser.add_argument("--output", "-out", required=True,
                       help="输出文件路径")
    parser.add_argument("--missing-output",
                       help="缺失基因输出文件路径")
    
    args = parser.parse_args()
    
    try:
        # 检查结果目录
        if not os.path.exists(args.results_dir):
            print(f"错误: 结果目录不存在: {args.results_dir}")
            sys.exit(1)
        
        # 合并结果
        print(f"正在合并结果目录: {args.results_dir}")
        merged_df = merge_results_from_directories(args.results_dir, args.output)
        
        print(f"成功合并 {len(merged_df)} 个探针")
        print(f"结果保存到: {args.output}")
        
        # 如果提供了基因列表，检查缺失的基因
        if args.gene_list:
            print(f"检查缺失基因: {args.gene_list}")
            
            # 加载原始基因列表
            gene_list = load_gene_list(args.gene_list)
            
            # 格式化基因名称
            formatted_genes = []
            for gene in gene_list:
                formatted_gene = format_gene_name(gene, args.organism)
                formatted_genes.append(formatted_gene)
            
            # 查找缺失的基因
            missing_genes = find_missing_genes(formatted_genes, merged_df)
            
            if missing_genes:
                print(f"发现 {len(missing_genes)} 个缺失基因:")
                for gene in missing_genes:
                    print(f"  - {gene}")
                
                # 保存缺失基因列表
                missing_output = args.missing_output or "missing_genes.txt"
                save_missing_genes(missing_genes, missing_output)
                print(f"缺失基因列表保存到: {missing_output}")
            else:
                print("所有基因都成功设计了探针！")
        
        # 显示统计信息
        print(f"\n统计信息:")
        print(f"  总探针数量: {len(merged_df)}")
        print(f"  涉及基因数量: {merged_df['gene_name'].nunique()}")
        print(f"  平均每个基因探针数: {len(merged_df) / merged_df['gene_name'].nunique():.1f}")
        
        if 'g_content' in merged_df.columns:
            print(f"  平均G含量: {merged_df['g_content'].mean():.2f}")
        if 'tm' in merged_df.columns:
            print(f"  平均熔解温度: {merged_df['tm'].mean():.1f} °C")
        
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
