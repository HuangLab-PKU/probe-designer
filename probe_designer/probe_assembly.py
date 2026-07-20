"""
Probe assembly module
Simplified probe assembly from binding sites with backbones.

Output Excel/JSON/FASTA follow the unified schema in
``probe_designer.io.probe_schema``: first three columns are
``order, probe_name, probe_sequence``, and ``probe_name`` follows
``{gene}_{position}_dRNA_{codebook}_{No.}`` (or ``..._iLock_...`` when
the iLock modification is applied). The intermediate ``binding_sites``
format keeps its short keys (``arm5``/``arm3``/``st``/``en``).
"""

import os
import json
from pathlib import Path
import pandas as pd
from typing import List, Dict, Any, Optional, Union

from .config import ProbeConfig
from .io.probe_schema import (
    CHEM_DRNA,
    CHEM_ILOCK,
    apply_final_column_order,
    make_padlock_name,
    resolve_codebook,
)


def assemble_plain_padlock(arm5: str, arm3: str, backbone: str) -> str:
    """Assemble a plain padlock probe — ``arm5 + backbone + arm3``.

    Used by mutation ``dRNA`` and ``cDNA`` chemistries, and by the TCR
    pipeline. iLock assembly is in
    :func:`probe_designer.ext.mutation.probe.assemble_ilock` because the
    invader linker geometry is chemistry-specific.

    Caller is responsible for casing — this function preserves whatever the
    arms / backbone come in as.
    """
    return arm5 + backbone + arm3


def load_binding_sites(path: Union[str, Path]) -> Dict[str, List[Dict[str, Any]]]:
    """Dispatch-on-extension loader for binding sites data.

    Accepts a path to either a pipeline-produced JSON
    (``{"gene_name": [{"arm_3prime": ..., "arm_5prime": ..., ...}, ...]}``)
    or a filtered-results XLSX with the historical column names
    (``gene_name``, ``arm3``, ``arm5``, ``st``, ``en``, ``g_content``,
    ``tm``, ``tm3``, ``tm5``, optional ``isoform_overlap_num``).

    Returns the canonical dict-of-lists shape used throughout the library.
    Raises ``FileNotFoundError``, ``ValueError`` on unsupported formats or
    missing XLSX columns.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Binding sites file not found: {path}")

    suffix = path.suffix.lower()
    if suffix == ".json":
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    if suffix in (".xlsx", ".xls"):
        df = pd.read_excel(path)
        required = {"gene_name", "arm3", "arm5", "st", "en", "g_content", "tm", "tm3", "tm5"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Excel file missing required columns: {sorted(missing)}"
            )

        out: Dict[str, List[Dict[str, Any]]] = {}
        for _, row in df.iterrows():
            gene = str(row["gene_name"])
            site = {
                "arm_3prime": str(row["arm3"]),
                "arm_5prime": str(row["arm5"]),
                "st": int(row["st"]),
                "en": int(row["en"]),
                "g_content": float(row["g_content"]),
                "tm": float(row["tm"]),
                "tm_3prime": float(row["tm3"]),
                "tm_5prime": float(row["tm5"]),
            }
            if "isoform_overlap_num" in df.columns and pd.notna(row.get("isoform_overlap_num")):
                site["isoform_overlap_num"] = row["isoform_overlap_num"]
            out.setdefault(gene, []).append(site)
        return out

    raise ValueError(
        f"Unsupported binding-sites file format: {suffix} (use .json or .xlsx)"
    )


class ProbeAssembler:
    """Assemble probes from binding sites with backbones."""
    
    def __init__(self, probe_config: ProbeConfig):
        self.probe_config = probe_config
        self.backbone_data = None
        if probe_config.backbone_file:
            self._load_backbone_data()
    
    def _load_backbone_data(self):
        """Load backbone data from Excel file."""
        if not os.path.exists(self.probe_config.backbone_file):
            raise FileNotFoundError(f"Backbone file not found: {self.probe_config.backbone_file}")
        
        # Read Excel file, expecting 'No.' and 'Sequence' columns
        df = pd.read_excel(self.probe_config.backbone_file)
        
        # Validate required columns
        if 'No.' not in df.columns or 'Sequence' not in df.columns:
            raise ValueError(f"Backbone file must contain 'No.' and 'Sequence' columns. Found: {df.columns.tolist()}")
        
        # Create mapping from No. to Sequence (convert No. to string for consistency)
        self.backbone_data = {str(no): seq for no, seq in zip(df['No.'], df['Sequence'])}
        
        # Validate no duplicate No.
        if len(self.backbone_data) != len(df):
            raise ValueError(f"Duplicate 'No.' values found in backbone file")
    
    def assemble_probes(
        self,
        binding_sites: Dict[str, List[Dict[str, Any]]],
        gene_info_file: str
    ) -> pd.DataFrame:
        """
        Assemble probes from binding sites using gene/bds info from Excel file.
        
        Args:
            binding_sites: Dict mapping gene_name to list of binding site dicts.
                          Each site must have: arm_3prime, arm_5prime, st, en, g_content, tm
            gene_info_file: Path to Excel file with 'gene_name' and 'No.' columns.
                          - Per-gene mode: one row per gene (same No. for all sites of that gene)
                          - Per-bds mode: one row per binding site (different No. for each site)
        
        Returns:
            DataFrame with columns: gene_name, probe_id, probe_seq, arm_5prime, arm_3prime,
                                   backbone, st, en, g_content, tm
        """
        if self.backbone_data is None:
            raise ValueError("Backbone data not loaded. Please provide backbone_file in ProbeConfig.")

        codebook = resolve_codebook(
            self.probe_config.codebook,
            Path(self.probe_config.backbone_file),
        )

        # Load gene_info from Excel file
        if not os.path.exists(gene_info_file):
            raise FileNotFoundError(f"Gene info file not found: {gene_info_file}")
        
        gene_info_df = pd.read_excel(gene_info_file)
        
        # Validate required columns
        if 'gene_name' not in gene_info_df.columns or 'No.' not in gene_info_df.columns:
            raise ValueError(f"Gene info file must contain 'gene_name' and 'No.' columns. Found: {gene_info_df.columns.tolist()}")
        
        # Normalize gene names: strip whitespace and convert to string
        gene_info_df['gene_name'] = gene_info_df['gene_name'].astype(str).str.strip()
        
        # Determine assembly mode: check if each gene has multiple rows
        gene_counts = gene_info_df['gene_name'].value_counts()
        is_per_bds_mode = (gene_counts > 1).any()
        
        # Build mapping from gene_name to No.(s) - use normalized names
        if is_per_bds_mode:
            # Per-bds mode: group by gene_name and create list of No.s
            gene_info = {}
            for gene_name, group in gene_info_df.groupby('gene_name'):
                gene_info[gene_name] = group['No.'].tolist()
        else:
            # Per-gene mode: one No. per gene
            gene_info = dict(zip(gene_info_df['gene_name'], gene_info_df['No.']))
        
        # Create a case-insensitive lookup map for better matching
        gene_info_lower = {k.strip().upper(): (k, v) for k, v in gene_info.items()}
        
        probe_records = []
        missing_genes = []
        
        for gene_name, sites in binding_sites.items():
            if not sites:
                continue
            
            # Normalize gene name for matching
            gene_name_normalized = gene_name.strip()
            gene_name_upper = gene_name_normalized.upper()
            
            # Try exact match first
            if gene_name_normalized in gene_info:
                probe_ids = gene_info[gene_name_normalized]
            # Try case-insensitive match
            elif gene_name_upper in gene_info_lower:
                original_name, probe_ids = gene_info_lower[gene_name_upper]
                print(f"警告: 基因名大小写不匹配，使用 '{original_name}' 匹配 '{gene_name}'")
            else:
                missing_genes.append(gene_name)
                continue
            
            probe_ids = gene_info[gene_name]
            
            # Check if it's per-gene or per-site mode
            if isinstance(probe_ids, (str, int, float)):
                # Per-gene mode: one No. for all sites
                probe_id = str(probe_ids)  # Convert to string for consistency
                if probe_id not in self.backbone_data:
                    raise ValueError(f"Probe ID '{probe_id}' not found in backbone file. Available IDs: {list(self.backbone_data.keys())}")

                backbone = self.backbone_data[probe_id]

                # Use same backbone for all sites of this gene
                for site in sites:
                    arm5 = site['arm_5prime'].lower()
                    arm3 = site['arm_3prime'].lower()
                    probe_seq = self._assemble_probe_sequence(arm5, arm3, backbone)
                    probe_records.append(self._build_record(
                        gene_name=gene_name, site=site, probe_id=probe_id,
                        probe_seq=probe_seq, backbone=backbone, codebook=codebook,
                    ))

            elif isinstance(probe_ids, list):
                # Per-site mode: one No. per binding site
                if len(probe_ids) != len(sites):
                    raise ValueError(f"Number of probe IDs ({len(probe_ids)}) does not match number of binding sites ({len(sites)}) for gene '{gene_name}'")

                for site, probe_id in zip(sites, probe_ids):
                    probe_id = str(probe_id)  # Convert to string for consistency
                    if probe_id not in self.backbone_data:
                        raise ValueError(f"Probe ID '{probe_id}' not found in backbone file. Available IDs: {list(self.backbone_data.keys())}")

                    backbone = self.backbone_data[probe_id]
                    probe_seq = self._assemble_probe_sequence(
                        site['arm_5prime'], site['arm_3prime'], backbone,
                    )
                    probe_records.append(self._build_record(
                        gene_name=gene_name, site=site, probe_id=probe_id,
                        probe_seq=probe_seq, backbone=backbone, codebook=codebook,
                    ))
            else:
                raise ValueError(f"Invalid gene_info format for gene '{gene_name}'. Expected str/int or list, got {type(probe_ids)}")
        
        # Report missing genes if any
        if missing_genes:
            available_genes = sorted(set(gene_info_df['gene_name'].tolist()))
            error_msg = (
                f"以下 {len(missing_genes)} 个基因在 gene_info 文件中未找到:\n"
                f"  缺失的基因: {', '.join(missing_genes)}\n"
                f"  gene_info 文件中的基因: {', '.join(available_genes[:20])}"
            )
            if len(available_genes) > 20:
                error_msg += f" ... (共 {len(available_genes)} 个)"
            raise ValueError(error_msg)
        
        df = pd.DataFrame(probe_records)
        return apply_final_column_order(df, kind="padlock")

    def _build_record(
        self, *, gene_name: str, site: Dict[str, Any], probe_id: str,
        probe_seq: str, backbone: str, codebook: str,
    ) -> Dict[str, Any]:
        """Schema-v2 record dict for one assembled probe."""
        position = int(site["st"])
        no = int(probe_id)
        record: Dict[str, Any] = {
            "order": pd.NA,
            "probe_name": make_padlock_name(
                name=gene_name, position=position, chemistry=CHEM_DRNA,
                codebook=codebook, no=no,
            ),
            "probe_sequence": probe_seq,
            "gene_name": gene_name, "clone_id": "",
            "chemistry": CHEM_DRNA,
            "position": position,
            "target_sequence": site.get("target_sequence", pd.NA),
            "probe_arm5": site["arm_5prime"],
            "probe_arm3": site["arm_3prime"],
            "probe_length": len(probe_seq),
            "No.": no, "codebook": codebook,
            "backbone_name": pd.NA,
            "backbone_sequence": backbone,
            # Schema-v2 keeps both metrics as distinct columns. The thermal
            # filter (Phase 0+) populates both; pre-Phase-0 binding-site dicts
            # only have g_content, so gc_content stays NaN for those.
            "gc_content": site.get("gc_content"),
            "g_content": site.get("g_content"),
            "tm": site.get("tm"),
            "tm_5prime": site.get("tm_5prime"),
            "tm_3prime": site.get("tm_3prime"),
            "free_energy": site.get("free_energy"),
        }
        if "isoform_overlap_num" in site:
            record["isoform_overlap_num"] = site["isoform_overlap_num"]
        return record

    def _assemble_probe_sequence(self, arm_5prime: str, arm_3prime: str, backbone: str) -> str:
        """
        Assemble probe sequence: arm_5prime + backbone + arm_3prime
        
        Args:
            arm_5prime: 5' arm sequence
            arm_3prime: 3' arm sequence
            backbone: Backbone sequence
        
        Returns:
            Complete probe sequence
        """
        return arm_5prime + backbone + arm_3prime
    
    def add_ilock_modification(self, probe_df: pd.DataFrame, position: str = "3prime") -> pd.DataFrame:
        """
        Add iLock modification to probes.

        Concatenates an iLock spacer to the existing probe and flips the
        ``chemistry`` column to ``"iLock"`` so the schema reflects the
        modification (no separate ``modification`` column is emitted —
        per Schema-v2 the chemistry value carries that signal).

        Args:
            probe_df: DataFrame with ``probe_sequence`` and ``chemistry`` columns.
            position: Where to add iLock — ``"3prime"``, ``"5prime"``, or ``"both"``.
        """
        ilock_seq = "CCTTCCTTCCTTCCTTCCTT"  # Typical iLock spacer

        probe_df = probe_df.copy()
        probe_df["chemistry"] = CHEM_ILOCK

        if position == "3prime":
            probe_df["probe_sequence"] = probe_df["probe_sequence"] + ilock_seq
        elif position == "5prime":
            probe_df["probe_sequence"] = ilock_seq + probe_df["probe_sequence"]
        elif position == "both":
            probe_df["probe_sequence"] = (
                ilock_seq + probe_df["probe_sequence"] + ilock_seq
            )
        else:
            raise ValueError(
                f"Invalid position '{position}'. Must be '3prime', '5prime', or 'both'"
            )

        # Probe length and probe_name need to track the modified sequence.
        probe_df["probe_length"] = probe_df["probe_sequence"].str.len()
        probe_df["probe_name"] = [
            make_padlock_name(
                name=str(r["gene_name"]), position=int(r["position"]),
                chemistry=CHEM_ILOCK, codebook=str(r["codebook"]), no=int(r["No."]),
            )
            for _, r in probe_df.iterrows()
        ]
        return probe_df

    def save_probes(self, probe_df: pd.DataFrame, output_dir: str):
        """Save probes to Excel, FASTA and JSON, plus stats."""
        os.makedirs(output_dir, exist_ok=True)

        excel_file = os.path.join(output_dir, "probes.xlsx")
        probe_df.to_excel(excel_file, index=False)

        fasta_file = os.path.join(output_dir, "probes.fasta")
        with open(fasta_file, "w") as f:
            for _, row in probe_df.iterrows():
                f.write(f">{row['probe_name']}\n{row['probe_sequence']}\n")

        json_file = os.path.join(output_dir, "probes.json")
        probe_dict = probe_df.to_dict("records")
        with open(json_file, "w", encoding="utf-8") as f:
            json.dump(probe_dict, f, indent=2, ensure_ascii=False, default=str)

        stats = {
            "total_probes": len(probe_df),
            "unique_genes": probe_df["gene_name"].nunique(),
            "unique_probe_nos": probe_df["No."].nunique(),
            "avg_probe_length": float(probe_df["probe_length"].mean()),
            "avg_gc_content": float(probe_df["gc_content"].mean()),
            "avg_tm": float(probe_df["tm"].mean()),
        }
        stats_file = os.path.join(output_dir, "probe_stats.json")
        with open(stats_file, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2, ensure_ascii=False)

        return probe_df
