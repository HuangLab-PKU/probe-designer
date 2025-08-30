from seqfold import dg, fold
import random
from Bio.Blast import NCBIXML


# barcode generator for PRISM


## random seq generator and filter
def random_seq_list(length=20, num=50):
    nucleotides = ["A", "T", "C", "G"]
    sequence = [
        "".join(random.choice(nucleotides) for _ in range(length)) for j in range(num)
    ]
    return sequence


def hum_dis(seq1, seq2):
    if len(seq1) != len(seq2):
        print("seq_not_match")
    else:
        cont = 0
        for char in range(len(seq1)):
            if seq1[char] != seq2[char]:
                cont += 1
        return cont


def dna_sec_struct(seq, temp=45):
    # Predict the minimum free energy
    mfe = dg(seq, temp=temp)
    # `fold` returns a list of `seqfold.Struct` from the minimum free energy structure
    structs = fold(seq, temp=temp)
    return mfe, structs


def thre_by_blast(file, thre=18):
    pos = []
    with open(file, "r") as blast_output:
        blast_records = NCBIXML.parse(blast_output)
        for blast_record in blast_records:
            save = True
            for alignment in blast_record.alignments:
                # Iterate over the high-scoring pairs (HSPs) in the alignment
                for hsp in alignment.hsps:
                    # print("HSP score:", hsp.score)
                    if hsp.score >= thre:
                        pos.append(False)
                        save = False
                        break
            if save:
                pos.append(True)
    return pos


## create seq lib and barcode lib based on grade of channels
import pandas as pd
import itertools


def create_seq_lib(
    seq_list,
    color_fraction={
        "Red": [_ / 4 for _ in range(5)],
        "Green": [_ / 2 for _ in range(3)],
        "Blue": [_ / 4 for _ in range(5)],
        "Yellow": [_ / 4 for _ in range(5)],
    },
):
    seq_lib = pd.DataFrame(columns=["seq", "color", "grade", "fraction"])
    color_list = []
    for color in color_fraction.keys():
        color_list += [color] * len(color_fraction[color])
    seq_lib["color"] = color_list

    if len(seq_lib["color"]) > len(seq_list):
        print("Seq Not Enough")
        return ValueError

    seq_lib["seq"] = seq_list[: len(seq_lib["color"])]

    fra = []
    grade = []
    for color in color_fraction.keys():
        fra += color_fraction[color]
        grade += [_ for _ in range(len(color_fraction[color]))]
    seq_lib["fraction"] = fra
    seq_lib["grade"] = grade

    return seq_lib


def create_barcode_lib(
    seq_lib,
    color_order=["Green", "Red", "Blue", "Yellow"],
    sum_num=5,
    sum_list=["Red", "Yellow", "Blue"],
):
    color_order_seq = [_ + "seq" for _ in color_order]
    barcode_lib = pd.DataFrame(columns=["barcode"] + color_order + color_order_seq)

    grade_list = {
        color: list(seq_lib[seq_lib.color == color].grade.unique())
        for color in color_order
    }
    barcode_lib[color_order] = list(itertools.product(*grade_list.values()))

    for i in range(len(barcode_lib)):
        grade = barcode_lib.loc[i, color_order]

        barcode_sub_list = [
            list(
                seq_lib[
                    (seq_lib.color == color_order[_]) & (seq_lib.grade == grade[_])
                ].seq
            )[0]
            for _ in range(len(color_order))
        ]
        barcode_lib.loc[i, color_order_seq] = barcode_sub_list
        barcode_lib.loc[i, "barcode"] = "".join(barcode_sub_list)

    if sum_num:
        barcode_lib["sum"] = barcode_lib[sum_list].sum(axis=1)
        barcode_lib = barcode_lib[(barcode_lib["sum"] == sum_num)]
        barcode_lib.set_index("barcode", inplace=True)

    return barcode_lib
