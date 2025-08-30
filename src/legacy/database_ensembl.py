from tqdm import tqdm
import requests

def ensembl_transcripts(gene_symbol="BRCA1", species='human'):
    """
    通过 Ensembl REST API 获取基因的转录本和外显子注释。
    """
    server = "https://rest.ensembl.org"
    lookup_ext = (
        f"/lookup/symbol/{species}/{gene_symbol}"
        f"?expand=1;transcripts=1"
    )
    headers = {"Content-Type": "application/json"}

    response = requests.get(server+lookup_ext, headers=headers)
    if response.ok:
        decoded = response.json()
        transcripts = decoded.get("Transcript")
        if transcripts: return transcripts
        else: print(f"No transcripts found for gene {gene_symbol}")
    else: print(f"Failed to fetch data: {response.status_code}, {response.text}")

# def ensembl_overlap_exons(gene_symbol="BRCA1", species="human", coord_system_version="GRCh38", tqdm_args={'position': 0, 'leave': True}):
#     """获取基因所有转录本及其外显子信息（优化版）"""
#     server = "https://rest.ensembl.org"
#     headers = {"Content-Type": "application/json"}
#     # 获取基因ID
#     lookup_ext = f"/lookup/symbol/{species}/{gene_symbol}?"
#     gene_data = requests.get(server+lookup_ext, headers=headers).json()
#     gene_id = gene_data["id"]

#     # 获取该基因区域所有转录本列表
#     tx_overlap_ext = f"/overlap/id/{gene_id}?feature=transcript"
#     transcripts = requests.get(server+tx_overlap_ext, headers=headers).json()

#     # 并行处理每个转录本
#     with tqdm(total=len(transcripts), desc=f"Fetching {gene_symbol} transcripts", **tqdm_args) as pbar:
#         for transcript in transcripts:
#             transcript_id = transcript['id']
#             transcript['exons'] = []          
#             try:
#                 # 获取精确外显子信息（单API调用）
#                 lookup_ext = f"/lookup/id/{transcript_id}?expand=1"
#                 tx_data = requests.get(server+lookup_ext, headers=headers).json()
                
#                 # 按链方向排序外显子
#                 exons = sorted(tx_data.get('Exon', []), key=lambda x: x['start'])
#                 transcript['exons'] = exons
            
#             except Exception as e: print(f"Error processing {transcript_id}: {str(e)}")
#             pbar.update(1)
#     return transcripts

def ensembl_fetch_sequence_region(chrom, start, end, species='human', coord_system_version="GRCh38"):
    """
    实际获取序列数据。
    """
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/{species}/{chrom}:{start}..{end}:1?"
    options = ";".join([
        "content-type=application/json", 
        f"coord_system_version={coord_system_version}"
    ])
    headers = {"Content-Type": "application/json"}
    
    response = requests.get(server + ext + options, headers=headers)
    if response.ok: return response.json()["seq"]
    else: response.raise_for_status()

def ensembl_fetch_sequence_batch(chromosome, intervals, species='human', coord_system_version="GRCh38"):
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

def ensembl_fetch_sequences(positions, coord_system_version = "GRCh38", gap=50):
    # Ensembl REST API URL for GRCh38 batch sequence fetching
    server = "https://rest.ensembl.org"
    sequences_info = []
    
    for position in tqdm(positions):
        # Adjust start and end for the extra 50 base pairs
        chromosome = position["chr"].replace("chr", "")
        st = int(position["start"])  # Ensure start is not less than 1
        en = int(position["end"])
        strand = int(position["strand"])
        ext = f"/sequence/region/human/{chromosome}:{st}..{en}:{strand}?"
        options = ";".join([
                "content-type=application/json",
                f"coord_system_version={coord_system_version}",
                f"expand_3prime={gap}",
                f"expand_5prime={gap}",
            ])
        response = requests.get(server + ext + options)

        if response.ok:
            decoder = response.json()
            decoder['gene'] = position['gene']
            decoder['chr'] = position['chr']
            decoder['start'] = position['start']
            decoder['end'] = position['end']
            sequences_info.append(decoder)
        else:
            response.raise_for_status()

    return sequences_info


# def ensembl_genename2seq(gene_symbol="BRCA1", species="human", seq_type="cds", 
#                          tqdm_args={'position': 0, 'leave': True}):
#     server = "https://rest.ensembl.org"
#     lookup_ext = f"/lookup/symbol/{species}/{gene_symbol}"
#     headers = {"Content-Type": "application/json"}
#     gene_id = requests.get(url=server+lookup_ext, headers=headers).json()["id"]

#     transcripts_ext = f"/overlap/id/{gene_id}?feature=transcript"
#     transcripts = requests.get(url=server+transcripts_ext, headers=headers).json()

#     # Get sequences for each transcript
#     transcripts_with_seq = []
#     with tqdm(total=len(transcripts), desc=f"{gene_symbol}", **tqdm_args) as pbar_task:
#         for transcript in transcripts:
#             seq_ext = f"/sequence/id/{transcript['id']}?type={seq_type}"
#             response = requests.get(server+seq_ext, headers=headers)
#             if response.ok:
#                 decoded = response.json()
#                 transcript["seq"] = decoded["seq"]
#                 if transcript: return transcripts_with_seq.append(transcript)
#                 else: print(f"No transcripts found for gene {gene_symbol}")
#             else: print(f"Failed to fetch data: {response.status_code}, {response.text}")
#             pbar_task.update(1)
#     return transcripts_with_seq

def ensembl_genename2seq(gene_symbol="BRCA1", species='human', seq_type='cds'):
    """
    通过 Ensembl REST API 获取基因的转录本和外显子注释。
    """
    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}

    # get transcript list for the gene
    transcripts = ensembl_transcripts(gene_symbol, species)

    # get sequences for each transcript
    for transcript in transcripts:
        seq_ext = f"/sequence/id/{transcript['id']}?type={seq_type}"
        seq_response = requests.get(server + seq_ext, headers=headers)
        if seq_response.ok: 
            transcript["seq"] = seq_response.json().get("seq", "")

    # filter out transcripts without sequence
    transcripts = [tx for tx in transcripts if 'seq' in tx]

    return transcripts


# def ensembl_geneid2seq(gene_id='ENSMUST00000128841.1', seq_type="cds", tqdm_args={'position': 0, 'leave': True}):
#     server = "https://rest.ensembl.org"
#     overlap_ext = (f"/overlap/id/{gene_id}?feature=transcript")
#     headers = {"Content-Type": "application/json"}
#     transcripts = requests.get(url=server+overlap_ext, headers=headers).json()
    
#     # Get sequences for each transcript
#     sequences = []
#     for transcript in tqdm(transcripts, desc=f"GeneID:\t{gene_id}"):
#         try:
#             seq_url = f"http://rest.ensembl.org/sequence/id/{transcript['id']}?type={seq_type};content-type=application/json"
#             seq_response = requests.get(seq_url).json()
#             transcript["seq"] = seq_response["seq"]
#             sequences.append(transcript)
#         except:
#             continue

#     return sequences

def ensembl_geneid2seq(transcript_id='ENSMUST00000128841.1', seq_type="cds", species='human', tqdm_args={'position': 0, 'leave': True}):
    """
    通过转录本ID获取对应基因的所有转录本序列
    步骤：
    1. 通过转录本ID查询基因信息
    2. 获取该基因所有转录本
    3. 获取每个转录本的指定类型序列
    
    参数：
    transcript_id: 输入转录本ID (isoform ID)
    seq_type: 序列类型 (cds, protein, genomic)
    species: 物种名称
    """
    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    
    # 第一步：通过转录本ID获取基因信息
    try:
        lookup_ext = f"/lookup/id/{transcript_id}?expand=0"
        tx_info = requests.get(server+lookup_ext, headers=headers).json()
        gene_id = tx_info['Parent']
        gene_symbol = tx_info.get('display_name', 'Unknown')
    except Exception as e:
        print(f"无法获取基因信息: {str(e)}")
        return []

    # 第二步：获取该基因所有转录本
    try:
        gene_lookup_ext = f"/lookup/symbol/{species}/{gene_symbol}?expand=1"
        gene_data = requests.get(server+gene_lookup_ext, headers=headers).json()
        transcripts = gene_data.get('Transcript', [])
    except Exception as e:
        print(f"无法获取转录本列表: {str(e)}")
        return []

    # 第三步：获取每个转录本序列
    sequences = []
    for tx in tqdm(transcripts, desc=f"Gene {gene_symbol}[{gene_id}]", **tqdm_args):
        try:
            seq_url = f"{server}/sequence/id/{tx['id']}?type={seq_type}"
            seq_data = requests.get(seq_url, headers=headers).json()
            sequences.append({
                "transcript_id": tx['id'],
                "gene_symbol": gene_symbol,
                "gene_id": gene_id,
                "seq_type": seq_type,
                "sequence": seq_data['seq'],
                "biotype": tx.get('biotype', '')
            })
        except Exception as e:
            print(f"跳过转录本 {tx['id']}: {str(e)}")
            continue
            
    return sequences

def get_strand_by_symbol(gene_symbol, species):
    """通过基因符号查询链信息"""
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/{species}/{gene_symbol}"
    headers = {"Content-Type": "application/json"}
    try:
        response = requests.get(server+ext, headers=headers)
        if response.ok: 
            return response.json().get('strand')
        else:
            print(f"Warning: {gene_symbol} not found (HTTP {response.status_code})")
            return None
    except Exception as e:
        print(f"Error for {gene_symbol}: {e}")
        return None

