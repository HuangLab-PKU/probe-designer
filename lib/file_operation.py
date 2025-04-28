# Env
import os
import pandas as pd

def merge(workdir=str, resultsdir="./results/", gene_name_file="immune.xlsx"):
    os.chdir(workdir)
    resultsdir = os.listdir("./results/")
    result = pd.DataFrame()
    for dir in resultsdir:
        try:
            df = pd.read_excel("./results/" + dir + "/probes_wanted.xlsx")
            result = pd.concat([result, df])
        except:
            continue

    result.to_excel("result.xlsx")

    gene_name = pd.read_excel(gene_name_file)
    gene_name_list = list(set(gene_name["gene_name"]))

    gene_out = pd.read_excel("result.xlsx")
    gene_out_name = list(set(gene_out["gene_name"]))

    gene_name_miss = [_ for _ in gene_name_list]
    for gene in gene_name_list:
        if gene in gene_out_name:
            gene_name_miss.remove(gene)

    with open("notfound.txt", "w") as f:
        f.write("\n".join(gene_name_miss))
