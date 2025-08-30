import numpy as np
from typing import Literal  # Python 3.8+
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

def get_color_list(n, cmap='viridis', format='hex', start=0.0, end=1.0):
    """修正后的颜色列表生成函数"""
    # 获取颜色映射对象
    cmap = plt.get_cmap(cmap)
    # 生成等间距位置
    positions = np.linspace(start, end, n)
    # 获取颜色数组（RGBA格式，0-1范围）
    colors = cmap(positions)
    # 格式转换
    if format == 'hex':
        return [mcolors.rgb2hex(color) for color in colors]  # 正确调用
    elif format == 'rgb':
        return [tuple(np.array(color[:3])*255).astype(int) for color in colors]
    elif format == 'rgba':
        return [tuple((np.array(color)*255).astype(int)) for color in colors]
    else:
        raise ValueError("Invalid format. Choose from ['hex', 'rgb', 'rgba']")

def get_style_config(is_unique, seq_type, style_type: Literal['unique', 'common', 'overlap', None] = None, overlap_num=0, global_overlap_num=1):
    colors = get_color_list(global_overlap_num, cmap='viridis', format='hex', start=0.0, end=1.0)
    """返回箭头样式配置字典"""
    style_config = {
        'exon': {
            'unique': {'color': 'red', 'label': 'Unique Exons', 'lw': 2.5},
            'common': {'color': 'orange', 'label': 'Common Exons', 'lw': 2},
            'overlap': {'color': colors[overlap_num], 'label': f'Overlapping Exons {overlap_num}', 'lw': 2.5},
        },
        'intron': {
            'unique': {'color': 'blue', 'label': 'Unique Introns', 'lw': 2.5},
            'common': {'color': 'lightblue', 'label': 'Common Introns', 'lw': 2}
        }
    }
    if style_type is None: style_type = 'unique' if is_unique else 'common'
    config = style_config[seq_type][style_type]
    # 修正箭头样式参数
    arrowstyle = "<->,head_length=0.2,head_width=0.15" if is_unique else "-"
    return {
        'arrowstyle': arrowstyle,
        'color': config['color'],
        'lw': config['lw'],
        'alpha': 0.9 if is_unique else 0.6,
        'label': config['label']
    }

def draw_strand_direction(ax, strand):
    """绘制链方向指示"""
    direction = 'right' if strand == 1 else 'left'
    pos_x = 0.1
    pos_y = 0.97
    arrow_length = 0.05
    if direction == 'right':
        arrow_end_x = pos_x + arrow_length
    else:
        arrow_end_x = pos_x - arrow_length
    ax.annotate('', xy=(arrow_end_x, pos_y), xytext=(pos_x, pos_y),
                xycoords='axes fraction', 
                arrowprops=dict(arrowstyle="->", color="red", lw=2))
    ax.text(0.02, 0.95, f"Transcription direction, Strand: {strand}",
            transform=ax.transAxes, color='red', fontsize=10)

# for unique regions in exons
def draw_transcript_unique(ax, tx, y_level):
    """绘制单个转录本的结构"""
    for seq_type, data_key, count_key in [('exon', 'exon_splices', 'exon_splices_count'), ('intron', 'introns', 'introns_count')]:
        for coord_pair in tx.get(data_key, []):
            # print(tx[count_key][coord_pair])
            style = get_style_config(tx[count_key][coord_pair]==1, seq_type)
            # 绘制箭头时只传递必要参数
            ax.annotate("", 
                        xy=(coord_pair[1], y_level), 
                        xytext=(coord_pair[0], y_level),
                        arrowprops={
                            'arrowstyle': style['arrowstyle'],
                            'color': style['color'],
                            'lw': style['lw'],
                            'alpha': style['alpha']
                        },
                        annotation_clip=False)
            # 添加标注
            if style['alpha'] == 0.9:
                text_offset = 0.4 if seq_type == 'exon' else -0.4
                ax.text(sum(coord_pair)/2, y_level + text_offset,
                        f"{coord_pair[0]}-{coord_pair[1]}", 
                        ha='center', va='center',
                        fontsize=9, color='darkred', weight='bold',
                        bbox=dict(boxstyle='round,pad=0.2', 
                                facecolor='#FFF8E1', 
                                alpha=0.8))
                
def plot_splice_structures(transcripts, figsize=None, ax=None):
    """主绘图函数"""
    if figsize is None:
        figsize = (2*len(transcripts), 10)
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
    
    # 设置坐标范围
    min_coord = min(tx['start'] for tx in transcripts) - 500
    max_coord = max(tx['end'] for tx in transcripts) + 500
    # 绘制转录方向指示
    draw_strand_direction(ax, transcripts[0]['strand'])
    # 绘制每个转录本
    for idx, tx in enumerate(reversed(transcripts)):
        y_level = idx * 2
        ax.set_yticks(range(0, len(transcripts)*2, 2))
        draw_transcript_unique(ax, tx, y_level)
        ax.set_yticklabels([t['external_name'] for t in transcripts[::-1]], fontsize=10)
    # 设置图例
    legend_elements = list({Line2D([0], [0], color=v['color'], lw=v['lw'], label=v['label']) 
                           for seq_type in ['exon', 'intron'] 
                           for v in [get_style_config(True, seq_type), get_style_config(False, seq_type)]})
    ax.legend(handles=legend_elements, loc='upper right')
    # 设置坐标轴
    ax.set_xlim(min_coord, max_coord)
    ax.set_ylim(-1, len(transcripts)*2)
    ax.set_xlabel("Genomic Position (chr5, GRCm39)")
    ax.set_title("Unique Splicing Patterns", pad=20)
    if ax is None:
        plt.grid(axis='x', linestyle=':')
        plt.tight_layout()
        plt.show()

# for overlap regions in exons
def draw_transcript_overlap(ax, tx, y_level):
    """绘制单个转录本的结构"""
    seq_type, data_key = 'exon', 'exon_splices_count'
    max_num = max(tx.get(data_key, {}).values(), default=1) + 1
    for coord_pair, num_overlap in tx.get(data_key, {}).items():
        style = get_style_config(False, seq_type, style_type='overlap', overlap_num=num_overlap, global_overlap_num=max_num)
        ax.annotate("", 
                    xy=(coord_pair[1], y_level), 
                    xytext=(coord_pair[0], y_level),
                    arrowprops={
                        'arrowstyle': style['arrowstyle'],
                        'color': style['color'],
                        'lw': style['lw'],
                        'alpha': style['alpha']
                    },
                    annotation_clip=False)
        # 添加标注
        if style['alpha'] == 0.9:
            text_offset = 0.4 if seq_type == 'exon' else -0.4
            ax.text(sum(coord_pair)/2, y_level + text_offset,
                    f"{coord_pair[0]}-{coord_pair[1]}", 
                    ha='center', va='center',
                    fontsize=9, color='darkred', weight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', 
                            facecolor='#FFF8E1', 
                            alpha=0.8))
            
def plot_overlap_structures(transcripts, figsize=None, ax=None):
    """主绘图函数"""
    if figsize is None:
        figsize = (2*len(transcripts), 10)
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
    # 设置坐标范围
    min_coord = min(tx['start'] for tx in transcripts) - 500
    max_coord = max(tx['end'] for tx in transcripts) + 500
    # 绘制转录方向指示
    draw_strand_direction(ax, transcripts[0]['strand'])
    # 绘制每个转录本
    for idx, tx in enumerate(reversed(transcripts)):
        y_level = idx * 2
        ax.set_yticks(range(0, len(transcripts)*2, 2))
        draw_transcript_overlap(ax, tx, y_level)
        ax.set_yticklabels([t['external_name'] for t in transcripts[::-1]], fontsize=10)
    # 设置图例
    max_overlap_num = max(tx.get('exon_splices_count', {}).values(), default=1)+1
    legend_elements = list({Line2D([0], [0], color=v['color'], lw=v['lw'], label=v['label']) 
                           for num_overlap in range(max_overlap_num)
                           for seq_type in ['exon'] 
                           for v in [get_style_config(False, seq_type, style_type='overlap', overlap_num=num_overlap, global_overlap_num=max_overlap_num)]})
    ax.legend(handles=legend_elements, loc='upper right')
    # 设置坐标轴
    ax.set_xlim(min_coord, max_coord)
    ax.set_ylim(-1, len(transcripts)*2)
    ax.set_xlabel("Genomic Position (chr5, GRCm39)")
    ax.set_title("Overlap Exons Patterns", pad=20)
    if ax is None:
        plt.grid(axis='x', linestyle=':')
        plt.tight_layout()
        plt.show()