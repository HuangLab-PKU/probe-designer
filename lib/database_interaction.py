from Bio import Entrez


# NCBI_database
def ncbi_get_GI(tmp, gene_name_list_file, gene_id_name_file):
    # Get gene id and other information from ncbi dataset(api)
    ## Generate gene_search_list from gene_name_list
    organism_of_interest = "Mus musculus"
    n_type_of_interest = "mRNA"
    with open(tmp + gene_name_list_file) as f:
        gene_name_list = f.read().splitlines()
    gene_search_list = [
        ", ".join([name, organism_of_interest, n_type_of_interest])
        for name in gene_name_list
    ]
    ## Get gene id list using Entrez.esearch
    id_list = []
    for gene_search in gene_search_list:
        Entrez.email = "1418767067@qq.com"
        handle = Entrez.esearch(db="nuccore", term=gene_search)
        record = Entrez.read(handle)
        handle.close()
        id_list += record["IdList"][:1]  # set number of search results to read
    with open(tmp + gene_id_name_file, "w") as f:
        f.write("\n".join(id_list))

def ncbi_get_genbank_from_GI(tmp, gene_id_name_file, genebank_file="3_gene_seq_in_file.gb"):
    ## Read id_list from existing file
    with open(tmp + gene_id_name_file, "r") as f:
        id_list = f.read().split("\n")

    # Get the genbank file of each gene by search for id list
    fetch_per_round = 3
    round = -(-len(id_list) // fetch_per_round)
    for i in range(round):
        id_list_per_round = id_list[i * fetch_per_round : (i + 1) * fetch_per_round]
        Entrez.email = "1418767067@qq.com"
        handle = Entrez.efetch(
            db="nuccore",
            strand=1,  # plus if strand=1
            id=id_list_per_round,
            rettype="gbwithparts",
            retmode="text",
        )
        seq_record = handle.read()
        handle.close()
        print(i + 1, "{:.2f} %".format((i + 1) / round * 100))
        with open(tmp + genebank_file, "a") as f:
            f.write(seq_record)


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

    big_seq = _fetch_sequence_region(chromosome.replace("chr",""), min_st, max_en, species, coord_system_version)
    seq_dict = {}
    for (st, en) in intervals:
        offset_s = st - min_st
        offset_e = en - min_st
        seq_dict[(st, en)] = big_seq[offset_s:offset_e]
    
    return seq_dict

def _fetch_sequence_region(chrom_clean, start, end, species='human', coord_system_version="GRCh38"):
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

# def ensembl_transcript_ex_in(gene="BRCA1", species="human", seq_type="cds", tqdm_args={'position': 0, 'leave': True}):
#     lookup_url = f"http://rest.ensembl.org/lookup/symbol/{species}/{gene}?content-type=application/json"
#     gene_id = requests.get(url=lookup_url).json()["id"]

#     # Get all transcripts for the gene
#     transcripts_url = f"http://rest.ensembl.org/overlap/id/{gene_id}?feature=transcript;content-type=application/json"
#     transcripts = requests.get(url=transcripts_url).json()

#     # Collect exons and introns for each transcript
#     # transcript_data = []
#     with tqdm(total=len(transcripts), desc=f"{gene}", **tqdm_args) as pbar_task:
#         for transcript in transcripts:
#             transcript_id = transcript['id']
#             try:
#                 seq_url = f"http://rest.ensembl.org/sequence/id/{transcript_id}?type={seq_type};content-type=application/json"
#                 seq_response = requests.get(seq_url).json()
#                 transcript["seq"] = seq_response["seq"]
            
#             except:
#                 pbar_task.update(1)
#                 continue

#             # Fetch exons for the transcript
#             exons_url = f"http://rest.ensembl.org/overlap/id/{transcript_id}?feature=exon;content-type=application/json"
#             exons = requests.get(url=exons_url).json()
#             # Sort exons by their start position
#             exons_sorted = sorted(exons, key=lambda x: x['start'])
#             # Deduce introns based on exons' positions
#             introns = []
#             for i in range(len(exons_sorted) - 1):
#                 introns.append({
#                     'intron_start': exons_sorted[i]['end'] + 1,
#                     'intron_end': exons_sorted[i + 1]['start'] - 1
#                 })

#             transcript['exons'] = exons_sorted,
#             transcript['introns'] = introns
#             pbar_task.update(1)

#     return transcripts
