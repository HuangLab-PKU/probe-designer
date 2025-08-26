def local_fetch_sequence(genome, seq_region, start, end):
    """从本地FASTA获取指定区域序列"""
    chrom = str(seq_region)  # 确保是字符串类型    
    # 验证染色体是否存在
    if chrom not in genome:
        raise KeyError(f"染色体 {chrom} 不存在于FASTA文件中，可用染色体: {list(genome.keys())[:10]}...")
    # 获取序列（自动处理负链）
    seq = genome[chrom][start-1:end]  # pyfaidx使用0-based坐标
    return str(seq)