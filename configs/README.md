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
  - `"brute_force"`: 暴力搜索，在整个序列上搜索所有可能的结合位点
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

```json
"filter": {
  "min_g_content": 0.3,
  "max_g_content": 0.7,
  "max_consecutive_g": 3,
  "min_tm": 60.0,
  "max_tm": 80.0,
  "max_tm_diff": 5.0,
  "min_free_energy": -5.0,
  "check_rna_structure": true,
  "max_alignments": 10,
  "require_specificity": true,
  "target_organisms": ["Mus musculus", "Homo sapiens"],
  "final_probes_per_gene": 3
}
```

#### min_g_content / max_g_content
- **含义**: G碱基含量的最小/最大值（0-1之间）
- **可能的值**: 0.0-1.0
- **推荐值**: `0.3-0.7`
- **默认值**: `0.3` / `0.7`

#### max_consecutive_g
- **含义**: 允许的最大连续G碱基数量
- **可能的值**: 正整数
- **推荐值**: `3-4`
- **默认值**: `3`

#### min_tm / max_tm
- **含义**: 熔解温度的最小/最大值（摄氏度）
- **可能的值**: 正数
- **推荐值**: `50-80°C`
- **默认值**: `60.0` / `80.0`

#### max_tm_diff
- **含义**: 3'和5'臂之间的最大熔解温度差异
- **可能的值**: 正数
- **推荐值**: `5-10°C`
- **默认值**: `5.0`

#### min_free_energy
- **含义**: RNA二级结构的最小自由能（kcal/mol）
- **可能的值**: 负数或0
- **推荐值**: `-10.0` 到 `0.0`
- **默认值**: `-5.0`

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
```json
{
  "filter": {
    "min_tm": 65.0,
    "max_tm": 75.0,
    "max_tm_diff": 3.0,
    "max_consecutive_g": 2
  }
}
```

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
