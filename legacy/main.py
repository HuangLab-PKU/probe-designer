# basci env
import os
import pandas as pd
import time
import json
from tqdm import tqdm

# data process of file from ncbi
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

# # get gene data from ncbi
# from Bio import Entrez

# # blast and xml file process
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# add package to sys var
# os.chdir(os.path.dirname(os.path.abspath(__file__)))
# sys.path.append("../lib")

# dir
workdir = './dataset/2024.3.27_mousebrain_HP_projection/'
os.makedirs(workdir, exist_ok=True)

current_time = time.localtime()
formatted_time = time.strftime("%Y%m%d_%H%M%S", current_time)

output = os.path.join(workdir, 'results', formatted_time)
pre_binding_dir = os.path.join(output, "pre_binding")
os.makedirs(output, exist_ok=True)

# basic variables
gene_name_list_tosearch = "gene_name_list_tosearch.txt"
pre_binding_file_suffix = "_pre_binding.fasta"
total_pre_binding_file_name = "_total.fasta"

# tmp file
pre_binding_num_file = "pre_binding_num.json"
blast_results_file = "blast_results.xml"


organism = 'mouse'
gene_info = pd.read_excel(os.path.join(workdir, "Gene list_mouse.xlsx"), sheet_name='Sheet1')
gene_list = list(gene_info['gene_name'].unique())


from lib.database_ensembl import ensembl_genename2seq
import time

max_trial = 3
sequences_of_all = dict()

# Dictionary to keep track of error messages for each gene
error_messages = {gene: [] for gene in gene_list}

with tqdm(total=len(gene_list), desc="total_process", position=0) as pbar_total:
    for gene in gene_list:
        sequences_of_all[gene] = {}
        retry_attempts = 0
        
        # Sub-task progress bar
        with tqdm(total=max_trial, desc=f"{gene}_trial", position=1, leave=False) as pbar_gene:
            while retry_attempts < max_trial:
                retry_attempts += 1
                try:
                    sequences_of_all[gene] = ensembl_genename2seq(gene=gene, species='mouse', seq_type='cds', tqdm_args={'position': 2,'leave': False})
                    break  # Successful retrieval, exit the retry loop
                
                except Exception as e:
                    # Log the error message and retrying info
                    tqdm.write(f"Error retrieving {gene}: {e}, retrying...", nolock=True)
                    time.sleep(1)  # Wait 1 second before retrying
                    
                    # Add error message to the dictionary
                    error_messages[gene].append(f"Attempt {retry_attempts}/{max_trial}: {e}")

                # Update the sub-task progress bar
                pbar_gene.update(1)
        
        # If all retry attempts failed, print the final error message
        if retry_attempts == max_trial:
            tqdm.write(f"Failed to retrieve sequences for {gene} after {max_trial} attempts.", nolock=True)
        
        # Update the main progress bar
        pbar_total.update(1)

# Print all error messages after processing all genes
for gene, messages in error_messages.items():
    for message in messages:
        print(message)