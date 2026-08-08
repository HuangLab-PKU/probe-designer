# 配置文件说明

本文档详细说明了`config_template.json`中各个配置选项的含义、可能的值和使用建议。

## 配置文件结构

配置文件采用JSON格式，包含以下主要部分：
- `database`: 数据库连接配置
- `search`: 结合位点搜索策略配置
- `filter`: 序列筛选条件配置
- `blast`: BLAST分析配置
- `output`: 输出文件配置
- `genome`: 基因组访问配置
- `species`: 物种配置

## 详细配置说明

### 1. Database 配置

```json
"database": {
  "database_type": "ensembl",
  "api_key": "your_ncbi_api_key_here",
  "email": "your_email@example.com"
}
```

#### database_type
- **含义**: 使用的数据库类型
- **可能的值**: 
  - `"ensembl"`: 使用Ensembl数据库（推荐）
  - `"ncbi"`: 使用NCBI数据库
- **默认值**: `"ensembl"`

#### api_key
- **含义**: NCBI API密钥（用于BLAST和数据库访问）
- **获取方式**: 在 https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ 申请
- **格式**: 字符串，如 `"your_ncbi_api_key_here"`

#### email
- **含义**: 用于数据库访问的邮箱地址
- **格式**: 有效的邮箱地址，如 `"your_email@example.com"`

### 2. Search 配置

```json
"search": {
  "search_strategy": "isoform_consensus",
  "binding_site_length": 40,
  "max_binding_sites": 30,
  "window_size": 50,
  "step_size": 1
}
```

#### search_strategy
- **含义**: 结合位点搜索策略
- **可能的值**:
  - `"single_sequence"`: 在一条输入序列上扫描所有合格位点（v0.2.0 前叫 `"brute_force"`）
  - `"exon_junction"`: 在外显子连接处搜索结合位点
  - `"isoform_specific"`: 搜索特定转录本的结合位点
  - `"isoform_consensus"`: 搜索覆盖多个转录本的共识结合位点（推荐）
- **默认值**: `"isoform_consensus"`

#### binding_site_length
- **含义**: 结合位点序列长度（碱基对）
- **可能的值**: 整数，通常为20-60
- **推荐值**: `40`（适合锁式探针）
- **默认值**: `40`

#### max_binding_sites
- **含义**: 每个基因的最大结合位点数量
- **可能的值**: 正整数
- **推荐值**: `30-100`
- **默认值**: `30`

#### window_size
- **含义**: 搜索窗口大小（用于某些搜索策略）
- **可能的值**: 正整数
- **默认值**: `50`

#### step_size
- **含义**: 搜索步长（用于暴力搜索）
- **可能的值**: 正整数
- **推荐值**: `1`（最精确）
- **默认值**: `1`

### 3. Filter 配置

> **Source of truth:** `probe_designer/config/FilterConfig` (defaults) and
> `configs/config_template.yaml` (annotated template). The values below are
> illustrative; when they disagree with the dataclass, the dataclass wins.

```yaml
filter:
  min_gc_content: 0.35       # Schema-v2 gate (GC%, not G-only)
  max_gc_content: 0.65
  max_consecutive_g: 5
  min_tm: 50.0               # arm-Tm target = lab_temp_c + tm_margin_c
  max_tm: 70.0
  max_tm_diff: 5.0           # two-arm balance tolerance
  enforce_tm_gate: false     # Tm rules are scoring references, not cutoffs
  enforce_tm_diff_gate: false
  min_free_energy: -10.0
  check_rna_structure: true
  max_alignments: 5
  require_specificity: true
  target_organisms: ["Mus musculus", "Homo sapiens"]
  final_probes_per_gene: 3
```

#### min_gc_content / max_gc_content
- **含义**: 每条臂的 G+C 含量下/上限（Schema-v2 实际生效的门槛）
- **推荐值**: `0.35-0.65`（2026-05-13 BZ23/TNBC 审计锚定）
- **默认值**: `0.35` / `0.65`
- 旧的 `min_g_content` / `max_g_content` 只数 G，已废弃，仅为向后兼容保留，不再生效。

#### max_consecutive_g
- **含义**: 允许的最大连续 G 碱基数量（A/T/C 同样有各自的门槛）
- **默认值**: `5`

#### min_tm / max_tm
- **含义**: 臂熔解温度的目标与上界（摄氏度），锚定到反应温度
  `lab_temp_c + tm_margin_c`，并在真实缓冲液（含 Mg²⁺ 与甲酰胺）下计算
- **默认值**: `50.0` / `70.0`
- **注意**: 默认是**打分参考**而非硬门槛 —— 离目标越远排名越低，但不会被丢弃。
  需要恢复硬性拒绝时把 `enforce_tm_gate` 设为 `true`。

#### max_tm_diff
- **含义**: 3' 与 5' 臂之间的熔解温度差异容忍度
- **推荐值**: `3-5°C`（Larsson 2010 / Krzywkowski 2017）
- **默认值**: `5.0`
- **注意**: 同样是打分参考；`enforce_tm_diff_gate` 才会让它变成硬门槛。在 232 条已下单
  探针上，5°C 的硬门槛会拒掉 17.7%，而作为打分参考则零代价。

#### min_free_energy
- **含义**: RNA 二级结构的最小自由能（kcal/mol）
- **可能的值**: 负数或0
- **默认值**: `-10.0`

#### check_rna_structure
- **含义**: 是否检查RNA二级结构
- **可能的值**: `true` / `false`
- **默认值**: `true`

#### max_alignments
- **含义**: BLAST比对的最大数量限制（已弃用，新版本不限制）
- **可能的值**: 正整数
- **默认值**: `10`

#### require_specificity
- **含义**: 是否要求特异性筛选
- **可能的值**: `true` / `false`
- **默认值**: `true`

#### target_organisms
- **含义**: 目标生物体列表（用于特异性筛选）
- **可能的值**: 字符串数组
- **示例**: `["Mus musculus", "Homo sapiens"]`
- **默认值**: `["Mus musculus", "Homo sapiens"]`

#### final_probes_per_gene
- **含义**: 每个基因的最终探针数量
- **可能的值**: 正整数
- **推荐值**: `3-5`
- **默认值**: `3`

### 4. BLAST 配置

```json
"blast": {
  "blast_type": "local",
  "database": "nt",
  "task": "megablast",
  "evalue": 1e-5,
  "hitlist_size": 50,
  "alignments": 500,
  "descriptions": 500,
  "megablast": true,
  "short_query": false,
  "filter": "none",
  "format_type": "XML",
  "service": "plain",
  "batch_size": 100,
  "concurrency": 3,
  "species": ["Mus musculus"]
}
```

#### blast_type
- **含义**: BLAST运行类型
- **可能的值**:
  - `"local"`: 本地BLAST（需要安装BLAST+）
  - `"online"`: 在线BLAST（使用NCBI服务器）
- **推荐值**: `"local"`（更快更稳定）
- **默认值**: `"local"`

#### database
- **含义**: BLAST数据库名称
- **可能的值**:
  - `"nt"`: 核苷酸数据库
  - `"refseq_rna"`: RefSeq RNA数据库（推荐）
  - `"nr"`: 非冗余蛋白质数据库
- **默认值**: `"nt"`

#### task
- **含义**: BLAST任务类型
- **可能的值**:
  - `"megablast"`: 快速核苷酸比对（推荐）
  - `"blastn"`: 标准核苷酸比对
  - `"blastp"`: 蛋白质比对
- **默认值**: `"megablast"`

#### evalue
- **含义**: E值阈值（显著性阈值）
- **可能的值**: 正数，通常为1e-3到1e-10
- **推荐值**: `1e-5`
- **默认值**: `1e-5`

#### hitlist_size
- **含义**: 返回的命中数量
- **可能的值**: 正整数
- **推荐值**: `50-100`
- **默认值**: `50`

#### alignments / descriptions
- **含义**: 显示的比对/描述数量
- **可能的值**: 正整数
- **推荐值**: `500`
- **默认值**: `500`

#### megablast
- **含义**: 是否使用megablast算法
- **可能的值**: `true` / `false`
- **默认值**: `true`

#### short_query
- **含义**: 是否为短查询调整参数
- **可能的值**: `true` / `false`
- **默认值**: `false`

#### filter
- **含义**: 序列过滤选项
- **可能的值**: `"none"`, `"L"`, `"R"`, `"F"`
- **默认值**: `"none"`

#### format_type
- **含义**: 输出格式
- **可能的值**: `"XML"`, `"Text"`
- **默认值**: `"XML"`

#### service
- **含义**: BLAST服务类型
- **可能的值**: `"plain"`, `"psi"`
- **默认值**: `"plain"`

#### batch_size
- **含义**: 在线BLAST的批次大小
- **可能的值**: 正整数
- **推荐值**: `100`
- **默认值**: `100`

#### concurrency
- **含义**: 并发请求数量
- **可能的值**: 正整数
- **推荐值**: `3-5`
- **默认值**: `3`

#### species
- **含义**: BLAST搜索的目标物种列表
- **可能的值**: 字符串数组
- **示例**: `["Mus musculus"]`, `["Mus musculus", "Homo sapiens"]`
- **默认值**: `["Mus musculus"]`
- **说明**: 如果未设置，将使用配置文件中的物种设置

### 5. Output 配置

```json
"output": {
  "output_dir": "results",
  "save_fasta": true,
  "save_json": true,
  "save_excel": true
}
```

#### output_dir
- **含义**: 输出目录路径
- **可能的值**: 字符串路径
- **默认值**: `"results"`

#### save_fasta
- **含义**: 是否保存FASTA格式文件
- **可能的值**: `true` / `false`
- **默认值**: `true`

#### save_json
- **含义**: 是否保存JSON格式文件
- **可能的值**: `true` / `false`
- **默认值**: `true`

#### save_excel
- **含义**: 是否保存Excel格式文件
- **可能的值**: `true` / `false`
- **默认值**: `true`

### 6. Genome 配置

```json
"genome": {
  "use_local_first": true
}
```

#### use_local_first
- **含义**: 是否优先使用本地基因组文件
- **可能的值**: `true` / `false`
- **默认值**: `true`

### 7. Species 配置

```json
"species": "mouse"
```

#### species
- **含义**: 目标物种名称
- **可能的值**: `"mouse"`, `"human"`, `"rat"`, `"zebrafish"`
- **默认值**: `"mouse"`

## 配置建议

### 小鼠基因组分析
```json
{
  "species": "mouse",
  "blast": {
    "species": ["Mus musculus"]
  },
  "filter": {
    "target_organisms": ["Mus musculus"]
  }
}
```

### 人类基因组分析
```json
{
  "species": "human",
  "blast": {
    "species": ["Homo sapiens"]
  },
  "filter": {
    "target_organisms": ["Homo sapiens"]
  }
}
```

### 多物种比较分析
```json
{
  "blast": {
    "species": ["Mus musculus", "Homo sapiens", "Rattus norvegicus"]
  },
  "filter": {
    "target_organisms": ["Mus musculus", "Homo sapiens", "Rattus norvegicus"]
  }
}
```

### 高特异性探针设计

收紧容忍度并把它们打开成硬门槛 —— 只有在你确实愿意为特异性丢掉候选时才这么做。

```yaml
filter:
  max_tm_diff: 3.0            # 更严格的双臂平衡
  max_consecutive_g: 4
  enforce_tm_gate: true       # 让 [min_tm, max_tm] 真正拒绝候选
  enforce_tm_diff_gate: true
```

若只是想提高杂交温度，改 `lab_temp_c` 即可 —— `min_tm` 会跟着
`lab_temp_c + tm_margin_c` 一起移动，不需要手动改阈值。

### 快速筛选模式
```json
{
  "search": {
    "max_binding_sites": 10
  },
  "filter": {
    "final_probes_per_gene": 1
  }
}
```
