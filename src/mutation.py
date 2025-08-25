import re
from tqdm import tqdm
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIXML
from copy import deepcopy

# update strand info
def extract_blast(sequences, blast_results):
    sequences_tmp = deepcopy(sequences)
    # read the id/plus-minus part/align_num
    with open(blast_results, "r") as blast_output:
        blast_records = NCBIXML.parse(blast_output)
        for num, blast_record in enumerate(blast_records):
            length = len(blast_record.alignments)
            sequences_tmp[num]["align_num"] = length
            sequences_tmp[num]['descrip'] = dict()
            for i in range(length):
                pm = blast_record.alignments[i].hsps[0].frame[1]
                sequences_tmp[num]["descrip"][i+1] = (blast_record.descriptions[i].title + f"| pm={pm}")
    return sequences_tmp


from lib.search_binding import seq_minus


def padlock_thre(Tm_left, Tm_right, left, right, Tm_dif_thre=10, Tm_thre_low=45, Tm_thre_high=65, bds_len_thre_low=10, bds_len_thre_high=30):
    flag = False
    if abs(Tm_left - Tm_right) < Tm_dif_thre:
        if Tm_left < Tm_thre_low: 
            left += 1
        elif Tm_left > Tm_thre_high:
            left -= 1
        elif Tm_right < Tm_thre_low: 
            right += 1
        elif Tm_right > Tm_thre_high:
            right -= 1
        else:
            flag = True

    elif Tm_left > Tm_right:
        left -= 1
        right += 1
        if left < bds_len_thre_low: flag = True
        if right > bds_len_thre_high: flag = True
        
    else:
        left += 1
        right -= 1
        if right < bds_len_thre_low: flag = True
        if left > bds_len_thre_high: flag = True
    return left, right, flag

def binding_searcher(seq_entry, gap=50, left_length=20, right_length=20, Tm_dif_thre=10, Tm_thre_low=45, Tm_thre_high=65, bds_len_thre_low=10, bds_len_thre_high=30):
    left = left_length
    right = right_length
    seq_entry = deepcopy(seq_entry)
    position_query = seq_entry['query']
    strand = seq_entry['strand']
    match = re.search(r':(\d+)\.\.(\d+):', position_query)
    query_st = int(match.group(1))
    query_en = int(match.group(2))
        
    seq = seq_entry["seq"]
    
    # if 'plp' not in seq_entry:
    #     seq_entry['plp'] = []
    seq_entry['plp'] = []

    plp_info = dict()
    while True:
        binding_left = seq[gap + 1 + query_en-query_st - left: gap + 1 + query_en-query_st]
        binding_right = seq[-gap : -gap + right]
        if strand == 1:
            Tm_left = mt.Tm_NN(binding_left, nn_table=mt.R_DNA_NN1)
            Tm_right = mt.Tm_NN(binding_right, nn_table=mt.R_DNA_NN1)
        elif strand == -1:
            Tm_left = mt.Tm_NN(seq_minus(binding_left), nn_table=mt.R_DNA_NN1)
            Tm_right = mt.Tm_NN(seq_minus(binding_right), nn_table=mt.R_DNA_NN1)

        left, right, flag = padlock_thre(Tm_left, Tm_right, left, right, 
                                         Tm_dif_thre=Tm_dif_thre, Tm_thre_low=Tm_thre_low, Tm_thre_high=Tm_thre_high, 
                                         bds_len_thre_low=bds_len_thre_low, bds_len_thre_high=bds_len_thre_high)
        if flag: break
    if strand == 1:
        # plp info
        plp_info["plp_left"] = seq_minus(binding_right)
        plp_info["plp_right"] = seq_minus(binding_left)
        plp_info["plp_bds"] = plp_info["plp_left"] + plp_info["plp_right"]
        plp_info["plp_Tm_left"] = Tm_right
        plp_info["plp_Tm_right"] = Tm_left
        plp_info["plp_Tm"] = mt.Tm_NN(plp_info["plp_bds"], nn_table=mt.R_DNA_NN1)
    elif strand == -1:
        plp_info["plp_left"] = binding_left
        plp_info["plp_right"] = binding_right
        plp_info["plp_bds"] = plp_info["plp_left"] + plp_info["plp_right"]
        plp_info["plp_Tm_left"] = Tm_left
        plp_info["plp_Tm_right"] = Tm_right
        plp_info["plp_Tm"] = mt.Tm_NN(plp_info["plp_bds"], nn_table=mt.R_DNA_NN1)
    seq_entry['plp'].append(plp_info)
    
    return seq_entry


import re
from Bio.SeqUtils import MeltingTemp as mt
from copy import deepcopy
from lib.search_binding import seq_minus

def mut_seq(seq, st, en, ref, alt):
    if ref == "-": mut_seq = seq[:en] + alt + seq[en:]
    elif alt == "-": mut_seq = seq[:st] + seq[en:]
    else: mut_seq = seq[:st] + alt + seq[en:]
    return mut_seq

def mut_binding(seq, ref, alt, st, en, left, right):
    if ref == "-":
        add = len(alt)
        if add < left: binding_left = seq[:en][add - left :] + alt
        else: binding_left = alt[-left:]
        binding_right = seq[en:][:right]

    elif alt == "-":
        binding_left = seq[:st][-left:]
        binding_right = seq[en:][:right]

    else:
        add = len(alt)
        if add < left: binding_left = seq[:st][add - left :] + alt
        else: binding_left = alt[-left:]
        binding_right = seq[en:][:right]
    
    return binding_left, binding_right

def padlock_thre(Tm_left, Tm_right, left, right, Tm_dif_thre=10, Tm_thre_low=45, Tm_thre_high=65, bds_len_thre_low=10, bds_len_thre_high=30):
    flag = False
    if abs(Tm_left - Tm_right) < Tm_dif_thre:
        if Tm_left < Tm_thre_low: 
            left += 1
        elif Tm_left > Tm_thre_high:
            left -= 1
        elif Tm_right < Tm_thre_low: 
            right += 1
        elif Tm_right > Tm_thre_high:
            right -= 1
        else:
            flag = True

    elif Tm_left > Tm_right:
        left -= 1
        right += 1
        if left < bds_len_thre_low: flag = True
        if right > bds_len_thre_high: flag = True
        
    elif Tm_left < Tm_right:
        left += 1
        right -= 1
        if right < bds_len_thre_low: flag = True
        if left > bds_len_thre_high: flag = True
    return left, right, flag

def perform_mutation(seq_entry, mutation, gap=50, left_length=20, right_length=20, Tm_dif_thre=10, Tm_thre_low=45, Tm_thre_high=65, bds_len_thre_low=10, bds_len_thre_high=30):
    # init for sequence
    strand = seq_entry['strand']
    seq = seq_entry['seq']
    position_query = seq_entry['query']
    match = re.search(r':(\d+)\.\.(\d+):', position_query)
    query_st = int(match.group(1))
    query_en = int(match.group(2))
    st = gap  # Ensure start is not less than 1
    en = int(query_en - query_st + gap + 1)

    # init for binding site
    left = left_length
    right = right_length

    # init for mutation
    mut_sequence = deepcopy(seq_entry)

    # perform mutation
    mut_info = deepcopy(mutation)
    ref, alt = mut_info['ref'], mut_info['alt']
    if strand == 1:
        mut_info['coding_seq'] = mut_seq(seq, st, en, ref, alt)
    elif strand == -1:
        # ref = seq_minus(ref)
        # alt = seq_minus(alt)
        mut_info['coding_seq'] = seq_minus(mut_seq(seq, st, en, ref, alt))
        
    while True:
        binding_left, binding_right = mut_binding(seq, ref, alt, st, en, left, right)
        if strand == 1:
            Tm_left = mt.Tm_NN(binding_left, nn_table=mt.R_DNA_NN1)
            Tm_right = mt.Tm_NN(binding_right, nn_table=mt.R_DNA_NN1)
        elif strand == -1:
            Tm_left = mt.Tm_NN(seq_minus(binding_left), nn_table=mt.R_DNA_NN1)
            Tm_right = mt.Tm_NN(seq_minus(binding_right), nn_table=mt.R_DNA_NN1)
            
        left, right, flag = padlock_thre(Tm_left, Tm_right, left, right, 
                                         Tm_dif_thre=Tm_dif_thre, Tm_thre_low=Tm_thre_low, Tm_thre_high=Tm_thre_high, 
                                         bds_len_thre_low=bds_len_thre_low, bds_len_thre_high=bds_len_thre_high)
        if flag: break

    if strand == 1:
        # reverse to padlock info
        mut_info['plp_bds'] = seq_minus(binding_left + binding_right)
        mut_info['plp_left'] = seq_minus(binding_right)
        mut_info['plp_right'] = seq_minus(binding_left)
        mut_info['plp_Tm'] = mt.Tm_NN(mut_info['plp_bds'], nn_table=mt.R_DNA_NN1)
        mut_info['plp_Tm_left'] = Tm_right
        mut_info['plp_Tm_right'] = Tm_left

    elif strand == -1:
        mut_info['plp_bds'] = binding_left + binding_right
        mut_info['plp_left'] = binding_left
        mut_info['plp_right'] = binding_right
        mut_info['plp_Tm'] = mt.Tm_NN(mut_info['plp_bds'], nn_table=mt.R_DNA_NN1)
        mut_info['plp_Tm_left'] = Tm_left
        mut_info['plp_Tm_right'] = Tm_right

    if 'mut' not in mut_sequence: mut_sequence['mut'] = [mut_info]
    else: mut_sequence['mut'].append(mut_info)

    return mut_sequence


import RNA

def ilock_params_cal(target_seq, link, len_5, len_3):
    target_5 = target_seq[link-len_5+1:link+1]
    target_3 = target_seq[link:link+len_3]
    ilock_bds_5 = seq_minus(target_5)
    ilock_bds_3 = seq_minus(target_3)
    ilock_bds = seq_minus(target_seq)
    return {
        'ilock_Tm': round(mt.Tm_NN(target_seq, nn_table=mt.R_DNA_NN1), 2), 
        "ilock_Tm5'": round(mt.Tm_NN(target_5, nn_table=mt.R_DNA_NN1), 2), 
        "ilock_Tm3'": round(mt.Tm_NN(target_3, nn_table=mt.R_DNA_NN1), 2), 
        'mfe': round(RNA.fold(target_seq)[1], 2),
        'ilock_bds': ilock_bds,
        "ilock_bds3'": ilock_bds_3,
        "ilock_bds5'": ilock_bds_5
        }


from lib.database_local import local_fetch_sequence

def local_pos_on_target_seq(pos5, pos3, seq, strand):
    if strand == 1:
        pos5_target = pos5
        pos3_target = pos3
        target_seq = seq
    elif strand == -1:
        pos5_target = len(seq) - pos3 - 1
        pos3_target = len(seq) - pos5 - 1
        target_seq = seq_minus(seq)
    return pos5_target, pos3_target, target_seq


import RNA
from Bio.SeqUtils import MeltingTemp as mt
from lib.search_binding import seq_minus, plp_wanted, plp_evaluation

def iter_plp(target_seq, pos5, pos3, plp_params):
    param_test = []
    ini_arm_len = plp_params.get('ini_arm_len', 20)
    min_arm_len = plp_params.get('min_arm_len', 15)
    max_arm_len = plp_params.get('max_arm_len', 25)
    for arm5_len in range(min_arm_len, max_arm_len+1):
        for arm3_len  in range(min_arm_len, max_arm_len+1):
            mfe = round(RNA.fold(target_seq[pos3: pos3+arm3_len] + target_seq[pos5-arm5_len: pos5])[1], 2)
            Tm5 = round(mt.Tm_NN(target_seq[pos5-arm5_len: pos5], nn_table=mt.R_DNA_NN1), 2)
            Tm3 = round(mt.Tm_NN(target_seq[pos3: pos3+arm3_len], nn_table=mt.R_DNA_NN1), 2)
            score = round(abs(Tm5-Tm3) + abs(arm5_len-ini_arm_len) + abs(arm3_len-ini_arm_len) - mfe, 2)
            param_test.append({
                'arm5_len': arm5_len, 'arm3_len': arm3_len, 
                'target_seq': target_seq[pos5-arm5_len: pos5] + " " + target_seq[pos3: pos3+arm3_len],
                "plp_bds5'": seq_minus(target_seq[pos5-arm5_len: pos5]), 
                "plp_bds3'": seq_minus(target_seq[pos3: pos3+arm3_len]),
                "plp_Tm5'": Tm5, "plp_Tm3'": Tm3, 'mfe': mfe, 'score': score})
    param_test = sorted(param_test, key=lambda x: x['score'])
    # param_return = [_ for _ in param_test if plp_wanted(_, plp_params)]
    for record in param_test: plp_evaluation(record, plp_params)
    param_return = [_ for _ in param_test]
    return param_return

def global_pos_plp(pos5, pos3, global_st, seq, strand):
    if strand == 1:
        pos5_global = global_st + pos5
        pos3_global = global_st + pos3
    elif strand == -1:
        pos5_global = global_st + len(seq) - pos3 - 1
        pos3_global = global_st + len(seq) - pos5 - 1
    return pos5_global, pos3_global

from lib.database_local import local_fetch_sequence

def plps_for_mut_item(mut_info_list, genome, plp_params, pad=50):
    mut_info_list = deepcopy(mut_info_list)    
    for item in tqdm(mut_info_list):
        global_mut_st, global_mut_en = item['Start'], item['End']
        global_st = max(0, global_mut_st - pad)
        global_en = global_mut_en + pad

        local_seq = local_fetch_sequence(genome, seq_region=item['Chr'].replace('chr', ''), start=global_st, end=global_en)
        local_mut_st, local_mut_en = global_mut_st - global_st, global_mut_en - global_st

        # DEL
        if item["Alt"] == "-": 
            item['mut_type'] = 'DEL'
            # ref
            item['plp_ref'] = []
            for pos5 in range(local_mut_st-1, local_mut_en+1):
                pos3 = pos5 + 1
                pos5_global, pos3_global = global_pos_plp(pos5, pos3, global_st, local_seq, item['strand'])
                ## find plps on target seq
                pos5_target, pos3_target, target_seq = local_pos_on_target_seq(pos5, pos3, local_seq, item['strand'])
                plp_info = iter_plp(target_seq, pos5_target+1, pos3_target, plp_params) # +1保证pos5能取到
                item['plp_ref'] += plp_info
                ## post process
                for plp in item['plp_ref']: plp["plp_pos5'"], plp["plp_pos3'"] = pos5_global, pos3_global
            item['plp_ref'] = sorted(item['plp_ref'], key=lambda x: x['score'])

            # mut
            local_mut_seq = local_seq[:local_mut_en] + item["Alt"] + local_seq[local_mut_en:]
            ## define pos of plp local and on the genome
            pos5, pos3 = local_mut_st - 1, local_mut_en + 1
            pos5_global, pos3_global = global_pos_plp(pos5, pos3, global_st, local_mut_seq, item['strand'])
            ## find plps on target seq
            pos5_target, pos3_target, target_seq = local_pos_on_target_seq(pos5, pos3-len(item['Alt']), local_mut_seq, item['strand'])
            plp_info = iter_plp(target_seq, pos5_target+1, pos3_target, plp_params) # +1保证pos5能取到
            item['plp_mut'] = plp_info
            ## post process
            for plp in item['plp_mut']: plp["plp_pos5'"], plp["plp_pos3'"] = pos5_global, pos3_global
            item['plp_mut'] = sorted(item['plp_mut'], key=lambda x: x['score'])


        # INS
        elif item["Ref"] == "-": 
            item['mut_type'] = 'INS'
            assert local_mut_st == local_mut_en

            # ref
            ## define pos of plp local and on the genome
            pos5, pos3 = local_mut_st, local_mut_en + 1
            pos5_global, pos3_global = global_pos_plp(pos5, pos3, global_st, local_seq, item['strand'])
            ## find plps on target seq
            pos5_target, pos3_target, target_seq = local_pos_on_target_seq(pos5, pos3, local_seq, item['strand'])
            plp_info = iter_plp(target_seq, pos5_target+1, pos3_target, plp_params) # +1保证pos5能取到
            item['plp_ref'] = plp_info
            ## post process
            for plp in item['plp_ref']: plp["plp_pos5'"], plp["plp_pos3'"] = pos5_global, pos3_global
            item['plp_ref'] = sorted(item['plp_ref'], key=lambda x: x['score'])

            # mut
            local_mut_seq = local_seq[:local_mut_st] + item["Alt"] + local_seq[local_mut_en:]
            pos5, pos3 = local_mut_st, local_mut_en + 1
            pos5_global, pos3_global = global_pos_plp(pos5, pos3, global_st, local_mut_seq, item['strand'])
            ## find plps on target seq
            item['plp_mut'] = []
            for ins in range(len(item["Alt"])+1):
                pos5 = local_mut_st + ins
                pos3 = local_mut_st + ins + 1
                pos5_target, pos3_target, target_seq = local_pos_on_target_seq(pos5, pos3, local_mut_seq, item['strand'])
                plp_info = iter_plp(target_seq, pos5_target+1, pos3_target, plp_params)
                for plp in plp_info: plp["plp_pos5'"], plp["plp_pos3'"] = f"{pos5_global}+INS{ins}", f"{pos3_global}+INS{ins}"
                item['plp_mut'] += plp_info
            item['plp_mut'] = sorted(item['plp_mut'], key=lambda x: x['score'])

        # REPLACE
        else: 
            if global_mut_st == global_mut_en: item['mut_type'] = 'SNP'
            else: item['mut_type'] = 'MNP'
            # ref
            item['plp_ref'] = []
            for pos5 in range(local_mut_st-1, local_mut_en+1):
                pos3 = pos5 + 1
                pos5_global, pos3_global = global_pos_plp(pos5, pos3, global_st, local_seq, item['strand'])
                ## find plps on target seq
                pos5_target, pos3_target, target_seq = local_pos_on_target_seq(pos5, pos3, local_seq, item['strand'])
                plp_info = iter_plp(target_seq, pos5_target+1, pos3_target, plp_params)
                for plp in plp_info: plp["plp_pos5'"], plp["plp_pos3'"] = pos5_global, pos3_global
                item['plp_ref'] += plp_info
            item['plp_ref'] = sorted(item['plp_ref'], key=lambda x: x['score'])

            # mut
            local_mut_seq = local_seq[:local_mut_st] + item["Alt"] + local_seq[local_mut_en+1:]
            pos5, pos3 = local_mut_st, local_mut_st + 1
            pos5_global, pos3_global = global_pos_plp(pos5, pos3, global_st, local_mut_seq, item['strand'])
            ## find plps on target seq
            item['plp_mut'] = []
            for rps in range(len(item['Alt'])+1):
                pos5 = local_mut_st + rps
                pos3 = local_mut_st + rps + 1
                pos5_target, pos3_target, target_seq = local_pos_on_target_seq(pos5, pos3, local_mut_seq, item['strand'])
                plp_info = iter_plp(target_seq, pos5_target+1, pos3_target, plp_params)
                for plp in plp_info: plp["plp_pos5'"], plp["plp_pos3'"] = f"{pos5_global}+RPS{rps}", f"{pos3_global}+RPS{rps}"
                item['plp_mut'] += plp_info
            item['plp_mut'] = sorted(item['plp_mut'], key=lambda x: x['score'])

    return mut_info_list