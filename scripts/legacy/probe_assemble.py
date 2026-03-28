#!/usr/bin/env python3
"""
探针拼接脚本

该脚本用于将结合位点与条形码序列拼接，生成最终的探针序列。
整合了notebook 3. probe_stitch.ipynb的逻辑。

作者: 重构版本
"""

import argparse
import pandas as pd
from pathlib import Path
from typing import Optional

# 添加src目录到路径
import sys
sys.path.append(str(Path(__file__).parent.parent))

# 延迟导入条形码工具，避免在仅使用BACKBONE模式时的硬依赖
def _lazy_load_barcode_data(panel: str, barcode_dir: Path) -> pd.DataFrame:
    from barcode import load_barcode_data  # type: ignore
    return load_barcode_data(panel, barcode_dir)


class ProbeStitcher:
    """探针拼接器"""
    
    def __init__(self, workdir: str, organism: str = "mouse"):
        """
        初始化探针拼接器
        
        Args:
            workdir: 工作目录路径
            organism: 物种名称
        """
        self.workdir = Path(workdir)
        self.organism = organism
        
    def load_binding_sites(self, binding_file: str = "gene_binding_site.xlsx") -> pd.DataFrame:
        """
        加载结合位点数据
        
        Args:
            binding_file: 结合位点文件路径
            
        Returns:
            结合位点DataFrame
        """
        binding_df = pd.read_excel(self.workdir / binding_file)
        
        # 根据物种格式化基因名称
        if self.organism == 'mouse':
            binding_df['gene_name'] = binding_df['gene_name'].str.capitalize()
        elif self.organism == 'human':
            binding_df['gene_name'] = binding_df['gene_name'].str.upper()
        
        print(f"加载了 {len(binding_df)} 个结合位点")
        return binding_df
    
    def load_gene_info(self, gene_info_file: str = "gene_info.csv") -> pd.DataFrame:
        """
        加载gene_info，并标准化关键字段以便映射到backbone。
        期望字段（宽松匹配，自动适配）：
        - 基因名列: gene | gene_name | Gene
        - 选择键列: id | probe_id | backbone_id | backbone_key
        - backbone序列列(可选): backbone | backbone_seq | sequence
        """
        path = self.workdir / gene_info_file
        if not path.exists():
            raise FileNotFoundError(f"未找到gene_info文件: {path}")
        if path.suffix.lower() in [".xlsx", ".xls"]:
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path)

        cols = {c.lower(): c for c in df.columns}
        # 基因列
        gene_col = None
        for k in ["gene", "gene_name", "Gene", "GENE"]:
            if k.lower() in cols:
                gene_col = cols[k.lower()]
                break
        if gene_col is None:
            raise ValueError("gene_info 中缺少基因列 (gene/gene_name)")
        df["gene"] = df[gene_col].astype(str)
        if self.organism == 'mouse':
            df["gene"] = df["gene"].str.capitalize()
        elif self.organism == 'human':
            df["gene"] = df["gene"].str.upper()

        # 选择键（用于选择backbone的id或key）
        select_key_col = None
        for k in ["id", "probe_id", "backbone_id", "backbone_key"]:
            if k in cols:
                select_key_col = cols[k]
                break
        if select_key_col is not None:
            df["backbone_key"] = df[select_key_col].astype(str)
        else:
            # 若没有显式key，则用基因名作为默认key
            df["backbone_key"] = df["gene"].astype(str)

        # 可能自带backbone序列
        for k in ["backbone", "backbone_seq", "sequence"]:
            if k in cols:
                df["backbone_seq"] = df[cols[k]]
                break
        return df[["gene", "backbone_key"] + (["backbone_seq"] if "backbone_seq" in df.columns else [])].drop_duplicates()

    def load_backbones(self, backbone_file: str) -> pd.DataFrame:
        """
        加载backbone库，支持CSV/XLSX。需要至少一个键列和一个序列列。
        键列优先: id | backbone_id | key | name
        序列列优先: backbone | sequence | seq | backbone_seq
        """
        path = Path(backbone_file)
        if not path.is_absolute():
            path = self.workdir / path
        if not path.exists():
            raise FileNotFoundError(f"未找到backbone文件: {path}")
        if path.suffix.lower() in [".xlsx", ".xls"]:
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path)
        cols = {c.lower(): c for c in df.columns}
        key_col = None
        for k in ["id", "backbone_id", "key", "name", "barcode", "label"]:
            if k in cols:
                key_col = cols[k]
                break
        if key_col is None:
            raise ValueError("backbone文件缺少键列(id/backbone_id/key/name)")
        seq_col = None
        for k in ["backbone", "sequence", "seq", "backbone_seq"]:
            if k in cols:
                seq_col = cols[k]
                break
        if seq_col is None:
            raise ValueError("backbone文件缺少序列列(backbone/sequence/seq)")
        out = df[[key_col, seq_col]].copy()
        out.columns = ["backbone_key", "backbone_seq"]
        out["backbone_key"] = out["backbone_key"].astype(str)
        return out.dropna(subset=["backbone_key", "backbone_seq"]).drop_duplicates()
    
    def load_barcode_data(self, panel: str = "PRISM", barcode_dir: str = None) -> pd.DataFrame:
        """
        加载条形码数据
        
        Args:
            panel: 面板类型 (PRISM/SPRINTseq)
            barcode_dir: 条形码文件目录
            
        Returns:
            条形码DataFrame
        """
        if barcode_dir is None:
            barcode_dir = self.workdir.parent
        
        barcode_df = _lazy_load_barcode_data(panel, barcode_dir)
        print(f"加载了 {len(barcode_df)} 个条形码")
        return barcode_df
    
    def stitch_prism_probes(self, binding_df: pd.DataFrame, barcode_df: pd.DataFrame) -> pd.DataFrame:
        """
        拼接PRISM探针
        
        Args:
            binding_df: 结合位点DataFrame
            barcode_df: 条形码DataFrame
            
        Returns:
            探针DataFrame
        """
        print("拼接PRISM探针...")
        
        probe_df = pd.DataFrame()
        
        # 获取PRISM条形码列表
        prism_list = [i for i in range(1, 31)]  # 可以根据需要调整
        
        for num, (prism, gene) in enumerate(zip(prism_list, binding_df["gene_name"].tolist())):
            if num >= len(binding_df):
                break
                
            binding = binding_df["bds"].iloc[num]
            binding_l = binding_df["binding_left"].iloc[num]
            binding_r = binding_df["binding_right"].iloc[num]
            
            # 获取条形码
            barcode_key = f"Prism_{prism}"
            if barcode_key not in barcode_df.index:
                print(f"警告: 未找到条形码 {barcode_key}")
                continue
                
            barcode = barcode_df.loc[barcode_key, "Barcode (82bp)"]
            
            # 拼接探针序列
            probe = binding_r.lower() + barcode.upper() + binding_l.lower()
            
            probe_info = pd.DataFrame({
                "PRISM": [f"PRISM_{prism}"],
                "gene": [f'{gene}'],
                "probe": [probe],
                "barcode": [barcode],
                "binding": [binding],
                "binding_left": [binding_l],
                "binding_right": [binding_r]
            })
            
            if len(probe_df) == 0:
                probe_df = probe_info
            else:
                probe_df = pd.concat([probe_df, probe_info])
        
        probe_df = probe_df.set_index('PRISM')
        return probe_df
    
    def stitch_sprintseq_probes(self, binding_df: pd.DataFrame, barcode_df: pd.DataFrame) -> pd.DataFrame:
        """
        拼接SPRINTseq探针
        
        Args:
            binding_df: 结合位点DataFrame
            barcode_df: 条形码DataFrame
            
        Returns:
            探针DataFrame
        """
        print("拼接SPRINTseq探针...")
        
        probe_df = pd.DataFrame()
        
        # SPRINTseq探针参数
        primer_l = 'TCCCTACACGACGCTCTTCCGATCT'
        primer_r = 'CATTCCTGCTGAACCGCTCTTCCGA'
        
        for num, gene in enumerate(binding_df["gene_name"].tolist()):
            binding = binding_df["bds"].iloc[num]
            binding_l = binding_df["binding_left"].iloc[num]
            binding_r = binding_df["binding_right"].iloc[num]
            
            # 获取条形码
            barcode_seq = barcode_df.iloc[num]["Barcode sequence"]
            barcode_70bp = primer_l + barcode_seq + primer_r + barcode_seq
            
            # 拼接探针序列
            probe = binding_r.lower() + barcode_70bp.upper() + binding_l.lower()
            
            probe_info = pd.DataFrame({
                "SPRINTseq": [f"SPRINTseq_{num+1}"],
                "gene": [f'{gene}'],
                "probe": [probe],
                "barcode": [barcode_70bp],
                "binding": [binding],
                "binding_left": [binding_l],
                "binding_right": [binding_r]
            })
            
            if len(probe_df) == 0:
                probe_df = probe_info
            else:
                probe_df = pd.concat([probe_df, probe_info])
        
        probe_df = probe_df.set_index('SPRINTseq')
        return probe_df
    
    def stitch_custom_probes(self, binding_df: pd.DataFrame, barcode_df: pd.DataFrame, 
                           panel: str, **custom_params) -> pd.DataFrame:
        """
        拼接自定义探针
        
        Args:
            binding_df: 结合位点DataFrame
            barcode_df: 条形码DataFrame
            panel: 面板类型
            **custom_params: 自定义参数
            
        Returns:
            探针DataFrame
        """
        print(f"拼接自定义 {panel} 探针...")
        
        # 这里可以根据需要实现其他类型的探针拼接
        # 暂时返回空DataFrame
        return pd.DataFrame()

    def stitch_backbone_probes(self, binding_df: pd.DataFrame, gene_info_df: pd.DataFrame, backbone_df: pd.DataFrame) -> pd.DataFrame:
        """
        按 arm5 + backbone + arm3 规则拼接通用探针。
        - 基于 gene_info 中的 backbone_key 选择对应的 backbone 序列
        - 若 gene_info 中已包含 backbone_seq，则优先使用；否则用 key 去 backbone_df 查找
        """
        print("拼接BACKBONE探针(arm5+backbone+arm3)...")
        # 标准化 binding_df 基因名
        work_df = binding_df.copy()
        if 'gene_name' not in work_df.columns:
            raise ValueError("binding表缺少gene_name列")
        work_df['gene'] = work_df['gene_name'].astype(str)
        if self.organism == 'mouse':
            work_df['gene'] = work_df['gene'].str.capitalize()
        elif self.organism == 'human':
            work_df['gene'] = work_df['gene'].str.upper()

        # 合并gene_info（带backbone_key/backbone_seq）
        merged = work_df.merge(gene_info_df, on='gene', how='left', suffixes=('', '_gi'))

        # 如 gene_info 未给出序列，使用 key 连接到 backbone 库
        need_lookup = 'backbone_seq' not in merged.columns
        if need_lookup or merged['backbone_seq'].isna().any():
            merged = merged.merge(backbone_df, on='backbone_key', how='left', suffixes=('', '_bb'))
            # 若存在两个backbone_seq列，优先gene_info提供的
            if 'backbone_seq_bb' in merged.columns:
                merged['backbone_seq'] = merged['backbone_seq'].combine_first(merged['backbone_seq_bb'])

        # 必要字段检查
        for col in ["binding_left", "binding_right", "backbone_seq"]:
            if col not in merged.columns:
                raise ValueError(f"拼接缺少必要列: {col}")

        # 构建探针
        probe_df = pd.DataFrame()
        for _, row in merged.iterrows():
            gene = row['gene']
            arm5 = str(row['binding_left'])
            arm3 = str(row['binding_right'])
            backbone_seq = str(row['backbone_seq'])
            if pd.isna(arm5) or pd.isna(arm3) or pd.isna(backbone_seq):
                continue
            probe_seq = arm5.lower() + backbone_seq.upper() + arm3.lower()
            rec = pd.DataFrame({
                "BACKBONE": [row.get('backbone_key', gene)],
                "gene": [gene],
                "probe": [probe_seq],
                "backbone": [backbone_seq],
                "binding_left": [arm5],
                "binding_right": [arm3]
            })
            if len(probe_df) == 0:
                probe_df = rec
            else:
                probe_df = pd.concat([probe_df, rec], ignore_index=True)
        if not probe_df.empty:
            probe_df = probe_df.set_index('BACKBONE')
        return probe_df
    
    def save_probes(self, probe_df: pd.DataFrame, panel: str):
        """
        保存探针序列
        
        Args:
            probe_df: 探针DataFrame
            panel: 面板类型
        """
        if probe_df.empty:
            print("警告: 没有生成任何探针")
            return
        
        # 保存Excel文件
        output_file = self.workdir / f"probes_{panel.lower()}.xlsx"
        probe_df.to_excel(output_file)
        print(f"保存探针到: {output_file}")
        
        # 保存FASTA文件
        fasta_file = self.workdir / f"probes_{panel.lower()}.fasta"
        with open(fasta_file, "w") as f:
            for idx, row in probe_df.iterrows():
                f.write(f">{idx}|{row['gene']}\n{row['probe']}\n")
        print(f"保存FASTA文件到: {fasta_file}")
        
        # 保存统计信息
        stats_file = self.workdir / f"probe_stats_{panel.lower()}.txt"
        with open(stats_file, "w") as f:
            f.write(f"面板类型: {panel}\n")
            f.write(f"探针总数: {len(probe_df)}\n")
            f.write(f"基因数: {probe_df['gene'].nunique()}\n")
            f.write(f"平均探针长度: {probe_df['probe'].str.len().mean():.1f} bp\n")
        print(f"保存统计信息到: {stats_file}")
    
    def run_stitch_pipeline(self, binding_file: str = "gene_binding_site.xlsx", 
                          panel: str = "PRISM", barcode_dir: str = None, 
                          backbone_file: Optional[str] = None, gene_info_file: str = "gene_info.csv",
                          **custom_params):
        """
        运行完整的探针拼接流程
        
        Args:
            binding_file: 结合位点文件
            panel: 面板类型
            barcode_dir: 条形码文件目录
            **custom_params: 自定义参数
        """
        print("开始运行探针拼接流程...")
        
        # 1. 加载结合位点
        binding_df = self.load_binding_sites(binding_file)
        
        if binding_df.empty:
            print("没有找到结合位点数据，流程结束")
            return
        
        # BACKBONE 模式无需条形码
        if panel.upper() == "BACKBONE":
            gene_info_df = self.load_gene_info(gene_info_file)
            backbone_df = self.load_backbones(backbone_file) if backbone_file else pd.DataFrame()
            probe_df = self.stitch_backbone_probes(binding_df, gene_info_df, backbone_df)
        else:
            # 2. 加载条形码数据
            barcode_df = self.load_barcode_data(panel, barcode_dir)
            if barcode_df.empty:
                print("没有找到条形码数据，流程结束")
                return
            # 3. 拼接探针
            if panel == "PRISM":
                probe_df = self.stitch_prism_probes(binding_df, barcode_df)
            elif panel == "SPRINTseq":
                probe_df = self.stitch_sprintseq_probes(binding_df, barcode_df)
            else:
                probe_df = self.stitch_custom_probes(binding_df, barcode_df, panel, **custom_params)
        
        # 4. 保存结果
        self.save_probes(probe_df, panel)
        
        print("探针拼接流程完成！")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="探针拼接工具")
    
    parser.add_argument("--workdir", required=True, help="工作目录路径")
    parser.add_argument("--binding_file", default="gene_binding_site.xlsx", help="结合位点文件")
    parser.add_argument("--organism", default="mouse", choices=["mouse", "human"], help="物种")
    parser.add_argument("--panel", default="PRISM", choices=["PRISM", "SPRINTseq", "BACKBONE"], help="面板类型/模式")
    parser.add_argument("--barcode_dir", help="条形码文件目录")
    parser.add_argument("--backbone_file", help="backbone文件路径(CSV/XLSX)")
    parser.add_argument("--gene_info_file", default="gene_info.csv", help="gene_info文件路径(CSV/XLSX)")
    
    args = parser.parse_args()
    
    # 创建探针拼接器
    stitcher = ProbeStitcher(args.workdir, args.organism)
    
    # 运行拼接流程
    stitcher.run_stitch_pipeline(
        binding_file=args.binding_file,
        panel=args.panel,
        barcode_dir=args.barcode_dir,
        backbone_file=args.backbone_file,
        gene_info_file=args.gene_info_file
    )


if __name__ == "__main__":
    main()

