from tqdm import tqdm
import requests

# Ensembl_database
def ensembl_name_to_seqs(gene="BRCA1", species="human", seq_type="cds", tqdm_args={'position': 0, 'leave': True}):
    lookup_url = f"http://rest.ensembl.org/lookup/symbol/{species}/{gene}?content-type=application/json"
    gene_id = requests.get(url=lookup_url).json()["id"]

    transcripts_url = f"http://rest.ensembl.org/overlap/id/{gene_id}?feature=transcript;content-type=application/json"
    transcripts = requests.get(url=transcripts_url).json()

    # Get sequences for each transcript
    sequences = []
    with tqdm(total=len(transcripts), desc=f"{gene}", **tqdm_args) as pbar_task:
        for transcript in transcripts:
            try:
                seq_url = f"http://rest.ensembl.org/sequence/id/{transcript['id']}?type={seq_type};content-type=application/json"
                seq_response = requests.get(seq_url).json()
                transcript["seq"] = seq_response["seq"]
                sequences.append(transcript)
                pbar_task.update(1)
            except:
                pbar_task.update(1)
                continue

    return sequences

def ensembl_id_to_seqs(gene="Gm16024", gene_id='ENSMUST00000128841.1', seq_type="cds"):
    transcripts_url = f"http://rest.ensembl.org/overlap/id/{gene_id}?feature=transcript;content-type=application/json"
    transcripts = requests.get(url=transcripts_url).json()

    # Get sequences for each transcript
    sequences = []
    for transcript in tqdm(transcripts, desc=f"Gene:\t{gene}"):
        try:
            seq_url = f"http://rest.ensembl.org/sequence/id/{transcript['id']}?type={seq_type};content-type=application/json"
            seq_response = requests.get(seq_url).json()
            transcript["seq"] = seq_response["seq"]
            sequences.append(transcript)
        except:
            continue

    return sequences

def ensembl_fetch_exons(gene_symbol, species='human', coord_system_version="GRCh38"):
    """
    通过 Ensembl REST API 获取基因的转录本和外显子注释。
    """
    server = "https://rest.ensembl.org"
    ext = (
        f"/lookup/symbol/{species}/{gene_symbol}"
        f"?expand=1;transcripts=1;coord_system_version={coord_system_version}"
    )
    headers = {"Content-Type": "application/json"}

    response = requests.get(server + ext, headers=headers)
    if response.ok:
        decoded = response.json()
        transcripts = decoded.get("Transcript")
        if transcripts:
            return transcripts
        else:
            print(f"No transcripts found for gene {gene_symbol}")
    else:
        print(f"Failed to fetch data: {response.status_code}, {response.text}")

def ensembl_fetch_sequence_region(chrom_clean, start, end, species='human', coord_system_version="GRCh38"):
    """
    实际获取序列数据。
    """
    server = "https://rest.ensembl.org"
    # 强制用正链(1)获取序列
    ext = f"/sequence/region/{species}/{chrom_clean}:{start}..{end}:1?"
    options = ";".join([
        "content-type=application/json",
        f"coord_system_version={coord_system_version}"
    ])
    headers = {"Content-Type": "application/json"}
    
    response = requests.get(server + ext + options, headers=headers)
    if response.ok:
        return response.json()["seq"]
    else:
        response.raise_for_status()

def ensembl_fetch_sequence_once(chromosome, intervals, species='human', coord_system_version="GRCh38"):
    """
    一次性获取多个区间的序列。
    
    参数:
        chromosome (str): 染色体名称，如 'chr11'
        intervals (list): [(start1, end1), (start2, end2), ...]
        
    返回:
        dict: {(start1, end1): sequence1, (start2, end2): sequence2, ...}
    """
    if not intervals:
        return {}
    
    min_st = min(iv[0] for iv in intervals)
    max_en = max(iv[1] for iv in intervals)

    big_seq = ensembl_fetch_sequence_region(chromosome.replace("chr",""), min_st, max_en, species, coord_system_version)
    seq_dict = {}
    for (st, en) in intervals:
        offset_s = st - min_st
        offset_e = en - min_st
        seq_dict[(st, en)] = big_seq[offset_s:offset_e]
    
    return seq_dict

def ensembl_transcript_exons(gene="BRCA1", species="human", seq_type="cds", tqdm_args={'position': 0, 'leave': True}):
    """获取基因所有转录本及其外显子信息（优化版）"""
    # 获取基因ID
    lookup_url = f"http://rest.ensembl.org/lookup/symbol/{species}/{gene}?content-type=application/json"
    gene_data = requests.get(lookup_url).json()
    gene_id = gene_data["id"]

    # 获取该基因所有转录本列表
    transcripts_url = f"http://rest.ensembl.org/overlap/id/{gene_id}?feature=transcript;content-type=application/json"
    transcripts = requests.get(transcripts_url).json()

    # 并行处理每个转录本
    with tqdm(total=len(transcripts), desc=f"Fetching {gene} transcripts", **tqdm_args) as pbar:
        for transcript in transcripts:
            transcript_id = transcript['id']
            # transcript['seq'] = None
            transcript['exons'] = []
            
            try:
                # 获取转录本序列
                # seq_url = f"http://rest.ensembl.org/sequence/id/{transcript_id}?type={seq_type};content-type=application/json"
                # seq_data = requests.get(seq_url).json()
                # transcript['seq'] = seq_data.get('seq')

                # 获取精确外显子信息（单API调用）
                lookup_url = f"http://rest.ensembl.org/lookup/id/{transcript_id}?expand=1;content-type=application/json"
                tx_data = requests.get(lookup_url).json()
                
                # 按链方向排序外显子
                exons = sorted(tx_data.get('Exon', []),
                              key=lambda x: x['start'], )
                            #   reverse=(tx_data['strand'] == -1))
                
                # 仅保留必要字段
                transcript['exons'] = exons
                
            except Exception as e:
                # 错误处理（跳过问题转录本）
                print(f"Error processing {transcript_id}: {str(e)}")
                pass
            
            pbar.update(1)

    # return [tx for tx in transcripts if tx['exons']]  # 过滤掉无外显子的转录本
    return transcripts  # 过滤掉无外显子的转录本