from copy import deepcopy
import re
from threading import local
from lib.search_binding import plp_search, seq_minus, position_search
from lib.database_local import local_fetch_sequence

# Version for unique plp_bds
def find_exon_counts(transcripts_raw):
    transcripts = deepcopy(transcripts_raw)
    # 收集所有exon及其出现次数
    exon_counts = {}
    # split exons
    split_points = sorted({point for tx in transcripts for exon in tx['exons'] for point in (exon['start'], exon['end'])})
    # 存储每个转录本的exons
    for transcript in transcripts:
        if not 'exons' in transcript: continue
        exons = transcript['exons']
        sorted_exons = sorted(exons, key=lambda x: x['start'])
        exon_splices = []
        # 遍历相邻的exon对生成exon
        for exon in sorted_exons:
            exon_start = exon['start']
            exon_end = exon['end']
            # use split points to split exon
            cur_splits = [exon_start] + [_ for _ in split_points if exon_start < _ < exon_end] + [exon_end]
            for j in range(len(cur_splits)-1):
                exon_splice = (cur_splits[j], cur_splits[j+1])
                exon_splices.append(exon_splice)
                exon_counts[exon_splice] = exon_counts.get(exon_splice, 0) + 1
        # save splices and junctions to transcript
        transcript['exon_splices'] = exon_splices
    for transcript in transcripts:
        transcript['exon_splices_count'] = {exon_splice: exon_counts[exon_splice] for exon_splice in transcript['exon_splices']}
    return transcripts

def find_intron_counts(transcripts_raw):
    transcripts = deepcopy(transcripts_raw)
    # 收集所有intron及其出现次数
    intron_counts = {}
    # 存储每个转录本的introns
    for transcript in transcripts:
        if not 'exons' in transcript: continue
        exons = transcript['exons']
        # 按rank升序排序
        sorted_exons = sorted(exons, key=lambda x: x['start'])
        introns = []
        # 遍历相邻的exon对生成intron
        for i in range(len(sorted_exons) - 1):
            current_exon = sorted_exons[i]
            next_exon = sorted_exons[i+1]
            intron = (current_exon['end'], next_exon['start'])
            introns.append(intron)
            # 更新全局计数
            intron_counts[intron] = intron_counts.get(intron, 0) + 1
        transcript['introns'] = introns
    for transcript in transcripts:
        transcript['introns_count'] = {intron: intron_counts[intron] for intron in transcript['introns']}
    return transcripts

def find_unique_introns(transcripts):
    # 确定每个转录本的独特intron
    unique_introns = []
    for transcript in transcripts:
        introns = transcript['introns']
        intron_counts = transcript['introns_count']
        unique = [j for j in introns if intron_counts[j] == 1]
        unique_introns.append({
            'transcript_id': transcript['id'], 
            'external_name': transcript['external_name'], 
            'intron': unique,
            })
    return unique_introns

def find_unique_exons(transcripts):
    # 确定每个转录本的独特exon
    unique_exons = []
    for transcript in transcripts:
        exons = transcript['exon_splices']
        exon_counts = transcript['exon_splices_count']
        unique = [j for j in exons if exon_counts[j] == 1]
        unique_exons.append({
            'transcript_id': transcript['id'], 
            'external_name': transcript['external_name'], 
            'exon_splices': unique,
            })
    return unique_exons

def find_overlap_exons(transcripts):
    transcripts_return = deepcopy(transcripts)
    # 确定每个转录本的独特exon
    overlap_exons = []
    for transcript in transcripts_return:
        overlap_exons.append({
            'transcript_id': transcript['id'], 
            'external_name': transcript['external_name'], 
            'exon_splices_count': transcript['exon_splices_count'],
            })
    return overlap_exons

def merge_adjacent_intervals(intervals):
    """合并相邻的连续区间"""
    if not intervals:
        return []
    # 按起始位置排序并转换类型便于修改
    sorted_intervals = sorted([list(interval) for interval in intervals], key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current[0] == last[1]:
            # 合并连续区间
            last[1] = current[1]
        else:
            merged.append(current)
    
    return [tuple(interval) for interval in merged]

def plp_bds_specific_regions(transcripts_raw, genome, pos_search_kwargs):
    transcripts = deepcopy(transcripts_raw)
    # 获取所有转录本的整体坐标范围
    global_start = min(t['start'] for t in transcripts)
    global_end = max(t['end'] for t in transcripts)
    # 计算全局序列中的坐标偏移
    def get_local_pos(pos): return pos - global_start

    # 获取整个基因座位的序列（假设所有转录本在同一染色体）
    sample_transcript = transcripts[0]
    seq = local_fetch_sequence(genome=genome, seq_region=sample_transcript['seq_region_name'], start=global_start, end=global_end)
    
    # 主处理逻辑
    for transcript in transcripts:
        # print(transcript['external_name'])
        transcript_id = transcript['id']
        regions = []
        transcript_probes = []

        # 处理exon_splices_count=1的区域
        candidate_splices = [
            splice for splice, count 
            in transcript['exon_splices_count'].items() 
            if count == 1
        ]
        merged_splices = merge_adjacent_intervals(candidate_splices)
        regions.extend(merged_splices)

        # 处理exon_splices_count=1的区域
        unique_regions = list(set(regions))
        final_regions = merge_adjacent_intervals(unique_regions)
        for region_start, region_end in final_regions:
            combined_seq = seq[get_local_pos(region_start): get_local_pos(region_end)]
            # 根据链方向处理序列拼接
            if transcript['strand'] == 1: rna_seq = combined_seq
            elif transcript['strand'] == -1: rna_seq = seq_minus(combined_seq)
            else: raise ValueError("Invalid strand value")
            # 设计探针
            pos_info = position_search(rna_seq, gene=transcript['external_name'], **pos_search_kwargs)
            for plp_probe in pos_info:
                if transcript['strand'] == 1: plp_probe['pos'] = global_start + plp_probe['pos']  
                elif transcript['strand'] == -1: global_start + len(rna_seq) - plp_probe['pos'] - pos_search_kwargs['BDS_len']
                else: raise ValueError("Invalid strand value")
            transcript_probes.extend(pos_info)

        # 处理introns_count=1的内含子两侧拼接
        for intron_idx, (intron, count) in enumerate(transcript['introns_count'].items()):
            if count == 1 and intron_idx < len(transcript['exons']) - 1:
                prev_exon = transcript['exons'][intron_idx]
                next_exon = transcript['exons'][intron_idx + 1]
                # 获取前exon末30bp在全局序列中的位置
                prev_start = max(prev_exon['start'], prev_exon['end'] - 30)
                prev_seq_start = get_local_pos(prev_start)
                prev_seq_end = get_local_pos(prev_exon['end'])
                # 获取后exon前30bp在全局序列中的位置
                next_seq_start = get_local_pos(next_exon['start'])
                next_seq_end = get_local_pos(min(next_exon['start'] + 30, next_exon['end']))
                # 根据链方向处理序列拼接
                combined_seq = seq[prev_seq_start:prev_seq_end] + seq[next_seq_start:next_seq_end]
                if transcript['strand'] == -1: rna_seq = seq_minus(combined_seq)
                else: rna_seq = combined_seq
                # 设计探针
                pos_info = position_search(rna_seq, gene=transcript['external_name'], **pos_search_kwargs)
                for plp_probe in pos_info:
                    # transform relative pos to global pos
                    if transcript['strand'] == 1: plp_probe['pos'] = global_start + plp_probe['pos']  
                    elif transcript['strand'] == -1: global_start + len(rna_seq) - plp_probe['pos'] - pos_search_kwargs['BDS_len']
                    else: raise ValueError("Invalid strand value")

                    # add overlap info:
                    plp_probe['isoform_overlap'] = 1
                transcript_probes.extend(pos_info)

        # 处理合并probes到transcript
        if 'plp_probes' in transcript: transcript['plp_probes'].extend(transcript_probes)
        else: transcript['plp_probes'] = transcript_probes
    
    return transcripts

# Version 2 for global target number
def find_exon_splice_counts(transcripts_raw):
    transcripts = deepcopy(transcripts_raw)
    # 收集所有exon及其出现次数
    exon_splice_counts = {}
    # split exons
    split_points = sorted({point for tx in transcripts for exon in tx['exons'] for point in (exon['start'], exon['end'])})
    # 存储每个转录本的exons
    for transcript in transcripts:
        transcript_name = transcript['external_name']
        if not 'exons' in transcript: continue
        exons = transcript['exons']
        sorted_exons = sorted(exons, key=lambda x: x['start'])
        # 遍历相邻的exon对生成exon
        transcript[f'exon_splices_1'] = []
        for exon in sorted_exons:
            exon_start = exon['start']
            exon_end = exon['end']
            # use split points to split exon
            cur_splits = [exon_start] + [_ for _ in split_points if exon_start < _ < exon_end] + [exon_end]
            for j in range(len(cur_splits)-1):
                exon_splice = (cur_splits[j], cur_splits[j+1])
                transcript[f'exon_splices_1'].append(exon_splice)
                exon_splice_counts.setdefault(exon_splice, []).append(transcript_name)
    for transcript in transcripts:
        transcript[f'exon_splices_1'] = {exon_splice: exon_splice_counts[exon_splice] for exon_splice in transcript[f'exon_splices_1']}
    return transcripts

def generate_min_sums(A):
    n = len(A)
    if n == 0:
        return []
    
    # 预处理前缀和数组
    prefix = [0] * (n + 1)
    for i in range(n):
        prefix[i + 1] = prefix[i] + A[i]
    
    # 初始化结果数组 B
    B = []
    for k in range(1, n + 1):  # k 代表窗口长度（从1到n）
        min_sum = float('inf')
        # 遍历所有长度为k的窗口
        for i in range(n - k + 1):
            # 计算窗口和：A[i] 到 A[i+k-1] 的和
            current_sum = prefix[i + k] - prefix[i]
            if current_sum < min_sum:
                min_sum = current_sum
        B.append(min_sum)
    
    return B

def exon_add_pos(start, length, splices):
    current_position = start  # 从给定的起始点开始
    remaining_length = length  # 还需覆盖的长度
    
    for i in range(len(splices)-1):
        start_point, end_point = splices[i]
        # 计算当前区间的有效长度（考虑起始点可能不在区间的开始处）
        if start_point <= current_position and current_position < end_point:
            remaining_in_interval = end_point - current_position
            if remaining_length < remaining_in_interval:
            # 如果剩余的长度可以在当前区间内完成，直接返回终点
                return current_position + remaining_length
            else:
            # 如果不能在当前区间完成，更新为下个区间的开始点
                remaining_length -= remaining_in_interval
                current_position = splices[i+1][0]
    return current_position + remaining_length

def exon_sub_pos(start, length, splices):
    current_position = start  # 从给定的起始点开始
    remaining_length = length  # 还需覆盖的长度
    
    # 反向遍历所有的splice区间
    for i in range(len(splices)-1, 0, -1):
        start_point, end_point = splices[i]
        # 如果当前起始点大于区间的结束点，更新为区间的结束点
        if current_position >= start_point and current_position < end_point:
            remaining_in_interval = current_position - start_point
            if remaining_length <= remaining_in_interval:
                # 如果剩余的长度可以在当前区间内完成，直接返回终点
                return current_position - remaining_length
            else:
                # 如果不能在当前区间完成，更新为下个区间的结束点
                remaining_length -= (remaining_in_interval + 1) 
                current_position = splices[i-1][1] - 1
    return current_position - remaining_length

def get_region(start, end, exon_splices):
    region = []
    for i in range(len(exon_splices)):
        exon_splice = exon_splices[i]
        if start >= exon_splice[0] and start < exon_splice[1]:
            region.append(exon_splice)
            if end >= exon_splice[0] and end < exon_splice[1]:
                break
            for j in range(i+1, len(exon_splices)):
                next_exon_splice = exon_splices[j]
                region.append(next_exon_splice)
                if end >= next_exon_splice[0] and end < next_exon_splice[1]:
                    break
    return tuple(region)

def get_seq(seq, global_start, start, end, region):
    def get_local_pos(pos): return pos - global_start
    # 获取目标序列
    target_seq = ''
    local_start = get_local_pos(start)
    local_end = get_local_pos(end)
    for exon_splice in region:
        exon_start = exon_splice[0]
        exon_end = exon_splice[1]
        # 获取在全局序列中的位置
        local_splice_start = get_local_pos(exon_start)
        local_splice_end = get_local_pos(exon_end)
        # 获取目标序列
        target_seq += seq[max(local_start, local_splice_start): min(local_end, local_splice_end)]
    return target_seq

from lib.search_binding import plp_params_cal

def plp_isoform(transcripts_raw, genome, plp_params):
    bds_len = plp_params.get('BDS_len')
    tlr_len = plp_params.get('junction_tolerence')

    transcripts = deepcopy(transcripts_raw)
    # 获取所有转录本的整体坐标范围
    global_start = min(t['start'] for t in transcripts)
    global_end = max(t['end'] for t in transcripts)

    # 获取整个基因区域的序列
    sample_transcript = transcripts[0]
    seq_region = sample_transcript['seq_region_name']
    seq = local_fetch_sequence(genome=genome, seq_region=seq_region, start=global_start, end=global_end)

    # calculate splice needed to define overlaps
    splice_def_needed = 1
    for transcript in transcripts:
        splice_length = [splice[1]-splice[0] for splice in transcript['exon_splices_1'].keys()]
        min_length = generate_min_sums(splice_length)
        for i in range(len(min_length)):
            if min_length[i] >= tlr_len:
                splice_def_needed = max(i + 2, splice_def_needed)
                break
    print(f"splice_def_needed: {splice_def_needed}")

    # generate splice combinations and calculate counts
    exon_splice_counts = {num: dict() for num in range(2, splice_def_needed+1)}
    for transcript in transcripts:
        exon_splices_1 = list(transcript['exon_splices_1'].keys())
        for num in range(2, splice_def_needed+1):
            transcript[f'exon_splices_{num}'] = []
        for num in range(2, splice_def_needed+1):
            for j in range(len(exon_splices_1)-num+1):
                exon_splice = tuple([exon_splices_1[j+k] for k in range(num)])
                transcript[f'exon_splices_{num}'].append(exon_splice)
                exon_splice_counts[num].setdefault(exon_splice, []).append(transcript['external_name'])

    # summary of all the transcripts
    for transcript in transcripts:
        for num in range(2, splice_def_needed+1):
            transcript[f'exon_splices_{num}'] = {exon_splice: exon_splice_counts[num][exon_splice] for exon_splice in transcript[f'exon_splices_{num}']}

    # walk on transcript to get probes
    for transcript in transcripts:
        transcript_start = transcript['start']
        transcript_end = transcript['end']
        exon_splices = list(transcript['exon_splices_1'].keys())

        transcript['padlock_probes'] = []
        start = transcript_start
        while True:
            start = exon_add_pos(start, 1, exon_splices)
            end = exon_add_pos(start, bds_len, exon_splices)
            if end >= transcript_end: break
            exon_splices_tmp = get_region(start, end, exon_splices)
            mid_point = exon_add_pos(start, bds_len // 2, exon_splices_tmp)

            # isoform judgement
            tlr_start = exon_sub_pos(mid_point, tlr_len, exon_splices_tmp)
            tlr_end = exon_add_pos(mid_point, tlr_len, exon_splices_tmp)
            tlr_region = get_region(tlr_start, tlr_end, exon_splices_tmp)
            # def the target isoform based on the tlr_regions
            exon_overlap_num = len(tlr_region)
            if exon_overlap_num == 1: 
                target_isoforms = transcript[f'exon_splices_{exon_overlap_num}'][tlr_region[0]]
            elif exon_overlap_num > 1: 
                target_isoforms = transcript[f'exon_splices_{exon_overlap_num}'][tlr_region]
            else:
                print(transcript['external_name'])
                print(start, tlr_start, mid_point, tlr_end, end)
                print(tlr_region)
                print(exon_splices_tmp)
                print(exon_splices)

            # get seq from start to end across exons
            region_seq = get_seq(seq, global_start, start, end, exon_splices_tmp)
            target_seq = region_seq if transcript['strand'] == 1 else seq_minus(region_seq)

            # get plp info from target_seq
            plp_info = {
                'start': start,
                'end': end,
                'exons': merge_adjacent_intervals(exon_splices_tmp),
                'target_seq': target_seq,
                'target_isoforms': target_isoforms,
                'exon_overlap_num': exon_overlap_num,
            }
            plp_info.update(plp_params_cal(target_seq, ))
            transcript['padlock_probes'].append(plp_info)

    return transcripts