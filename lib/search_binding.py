import os
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# import RNA
# from seq_filter import pre_blast_select, pre_box_select
from Bio.SeqUtils import MeltingTemp as mt


def pre_filter(gene_name, gene_name_list, mol_type, organism):
    select = True
    # check gene name
    if gene_name not in gene_name_list:
        select = False
    # get molecule_type
    if select:
        if mol_type != "mRNA":
            select = False
    # # check organism
    # if select:
    #     if organism != "Homo sapiens":
    #         select = False
    return select


def pre_blast_filter(seq):
    select = True
    # non consective 5 base
    if "GGGGG" in seq:
        select = False
    # check G 40%-70%
    if select == True:
        G_per = seq.count("G") / len(seq)
        if G_per < 0.3 or G_per > 0.7:
            select = False
    return select


def seq_minus(seq):
    translib = {"A": "T", "T": "A", "C": "G", "G": "C", "-":"-"}
    return "".join(list(reversed([translib[i] for i in seq])))


def gb_extract(record, gene_name='Nan', CDS=True):
    # get information of gene
    translib = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    coding_sequence = ""
    gene_id = record.id
    mol_type = record.annotations["molecule_type"]
    organism = record.annotations["organism"]

    if record.features:
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", [gene_name])[0]
                coding_sequence = feature.location.extract(record).seq

    if CDS:
        if coding_sequence == "":
            print(f"No CDS found in genbank file, please use another file or select on all sequence.")
        seq = coding_sequence
        # seq_minus = [translib[i] for i in str(coding_sequence)]
        # seq = "".join(list(reversed(seq_minus)))

    else:
        seq = record.seq
    #     seq_minus = [translib[i] for i in str(record.seq)]
    #     seq = "".join(list(reversed(seq_minus)))
    gb_info = {
        'gene_id': gene_id,
        'gene_name': gene_name,
        'mol_type': mol_type,
        'organism': organism,
        'seq': seq,
    }
    return gb_info


def site_searcher(
    tmp_output_pd,
    genbank_file=str,
    tmp="./tmp/",
    BDS_file_out_dir="./pre_binding_dir",
    pre_binding_file_suffix="_pre_binding.fasta",
    total_pre_binding_file_name="_total_pre_binding.fasta",
    pre_binding_num_file="4_pre_binding_num.json",
    gene_name_list=[],
    BDS_len=40,
    max_num=30,
):
    os.makedirs(BDS_file_out_dir, exist_ok=True)
    tmp_max_num = max_num * 2
    pre_binding_num = {}
    pos, length = 0, 0
    translib = {"A": "T", "T": "A", "C": "G", "G": "C"}

    with open(BDS_file_out_dir + total_pre_binding_file_name, "w") as handle:
        handle.write("")

    for gene_seq_in in SeqIO.parse(tmp + genbank_file, "genbank"):
        # get information of gene
        id = gene_seq_in.id  # get seq id
        # get gene_name
        if gene_seq_in.features:
            for feature in gene_seq_in.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ["NAN"])[0]

        # get molecule_type
        mol_type = gene_seq_in.annotations["molecule_type"]

        # get organism
        organism = gene_seq_in.annotations["organism"]

        # select by gene_name, mol_type and organism
        if not pre_filter(gene_name, gene_name_list, mol_type, organism): continue

        # get minus seq
        seq = str(gene_seq_in.seq)
        try: seq_minus = "".join(list(reversed([translib[i] for i in seq])))
        except: continue  # pass the loop if seq is not a mRNA seq

        # set start point and pre_binding_num
        length = len(seq_minus)
        pre_binding_num_tmp = min(length // BDS_len, tmp_max_num)
        st = length // 2 - pre_binding_num_tmp * BDS_len // 2

        valid_num = 0
        # generate and write pre_binding for each gene in a fasta file
        record_list = []
        file_out = f"{BDS_file_out_dir}{id}{pre_binding_file_suffix}"

        for i in range(pre_binding_num_tmp):
            pre_binding_tmp = seq_minus[st + i * BDS_len : st + i * BDS_len + BDS_len]

            if not pre_blast_filter(pre_binding_tmp): continue
            valid_num += 1

            record_list.append(
                SeqRecord(Seq(pre_binding_tmp),
                          id="pre_binding" + str(i),
                          description="|".join([id, gene_name, organism, mol_type])))

            # add information about binding sites to tmp_output_pd
            tmp_output_pd.loc[i + pos, ["accession", "gene_name", "mol_type", "organism", "binding"]
            ] = [id, gene_name, mol_type, organism, pre_binding_tmp]

            if valid_num >= max_num:
                break

        pos += valid_num

        # write pre_binding to files
        with open(file_out, "w") as f:
            for new_record in record_list:
                SeqIO.write(new_record, f, "fasta")
        with open(BDS_file_out_dir + total_pre_binding_file_name, "a") as handle:
            for new_record in record_list:
                SeqIO.write(new_record, handle, "fasta")

        # record the num of pre_binding for each gene
        pre_binding_num[id] = [gene_name, valid_num]

    with open(tmp + pre_binding_num_file, "w") as f:
        json.dump(pre_binding_num, f)

    return tmp_output_pd


import random
from tqdm import tqdm


def random_no_overlap_seqs(input_string, length, num_substrings):
    if length > len(input_string) or num_substrings * length > len(input_string):
        # if length > len(input_string):
        raise ValueError("Invalid input parameters.")

    available_positions = list(range(len(input_string) - length + 1))
    random.shuffle(available_positions)

    substrings = []
    positions = set()  # To keep track of selected positions

    for _ in range(num_substrings):
        if not available_positions:
            break  # Stop if there are no more non-overlapping positions

        # Try to find a non-overlapping position
        selected_position = available_positions.pop()
        while selected_position in positions:
            if not available_positions:
                break
            selected_position = available_positions.pop()

        if selected_position not in positions:
            positions.add(selected_position)
            end = selected_position + length
            selected_substring = input_string[selected_position:end]
            substrings.append(selected_substring)

    return substrings


def optimize_subsequence(
    positions, length, min_gap, better_gap=80,
    gene="", warn=True,):
    
    if len(positions) == 0:
        if warn: print(f"Gene {gene}: No valid positions, please loosen your threshold conditions.")
        return []

    positions.sort()

    def is_valid(min_difference, length, min_gap):
        count = 1
        current_min = positions[0]

        for i in range(1, len(positions)):
            if positions[i] - current_min >= min_difference:
                count += 1
                current_min = positions[i]

        return count >= length and min_difference > min_gap

    left, right = 0, positions[-1] - positions[0]
    result = []
    while left <= right:
        mid = (left + right) // 2
        if is_valid(mid, length, min_gap):
            result = [positions[0]]
            current_min = positions[0]
            for i in range(1, len(positions)):
                if positions[i] - current_min >= mid:
                    result.append(positions[i])
                    current_min = positions[i]
            left = mid + 1
        else:
            right = mid - 1

    if result == []:
        if warn: print(f"Gene {gene}: Not enough pos for {length} binding sites.")
        result = positions

    if mid < better_gap and warn:
        print(f"Gene {gene}: condition too harsh, loose to get better results")
        print(result)

    return result


import RNA

def position_search(sequence,  
    BDS_len=40, BDS_num=3, min_gap=0, better_gap=40, pin_gap=0.1, 
    G_min=0.3, G_max=0.7, G_consecutive=5, Tm_low=50, Tm_high=65, Tm_dif_thre=5, mfe_thre=-10, 
    gene="", verbose=True, verbose_pos=0, leave=True, warn=True):
    
    if len(sequence) < BDS_len:
        if warn: print(f"Gene {gene}: sequence too short, please use longer sequence.")
        return []
    
    seq_gap = int(len(sequence) * pin_gap)
    position = [_ for _ in range(seq_gap, len(sequence) - seq_gap - BDS_len)]
    plp_pos_candidate = []

    for pos in tqdm(position, desc=f"{gene}", disable=not verbose, position=verbose_pos, leave=leave):
        target_seq = sequence[pos : pos + BDS_len]
        # non consective base
        if "G" * G_consecutive in target_seq: continue
        # check G 40%-70%,
        G_pct = target_seq.count("G") / len(target_seq)
        if G_pct < G_min or G_pct > G_max: continue
        # check Tm
        plp_Tm_3 = mt.Tm_NN(target_seq[BDS_len//2:], nn_table=mt.R_DNA_NN1)
        plp_Tm_5 = mt.Tm_NN(target_seq[:BDS_len//2], nn_table=mt.R_DNA_NN1)
        plp_Tm = mt.Tm_NN(target_seq, nn_table=mt.R_DNA_NN1)
        if plp_Tm_3 > Tm_high or plp_Tm_3 < Tm_low or plp_Tm_5 > Tm_high or plp_Tm_5 < Tm_low: continue # Tm too high or too low
        if abs(plp_Tm_5 - plp_Tm_3) > Tm_dif_thre: continue # Tm difference too great
        # check 2nd structure
        _, mfe = RNA.fold(target_seq)
        if mfe < mfe_thre: continue

        plp_pos_candidate.append({
            'pos': pos, 
            'plp_Tm': round(plp_Tm,2), "plp_Tm5'": round(plp_Tm_5,2), "plp_Tm3'": round(plp_Tm_3,2), 
            'mfe': round(mfe,2)})

    pos_best = optimize_subsequence(
        [_['pos'] for _ in plp_pos_candidate], 
        BDS_num, min_gap=min_gap, better_gap=better_gap,
        gene=gene, warn=warn,)
    
    pos_info = []
    for record in plp_pos_candidate:
        if record['pos'] in pos_best:
            pos = record['pos']
            record['plp_bds'] = seq_minus(sequence[pos : pos + BDS_len])
            record["plp_bds3'"] = record['plp_bds'][: BDS_len // 2]
            record["plp_bds5'"] = record['plp_bds'][BDS_len // 2 :]
            pos_info.append(record)

    return pos_info
