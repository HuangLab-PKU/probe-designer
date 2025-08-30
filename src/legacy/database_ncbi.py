from Bio import Entrez
Entrez.email = "1418767067@qq.com"
Entrez.api_key = '010eacb785458478918b0cb14bea9f9df609'

# NCBI_database
def ncbi_get_gi(tmp, gene_name_list_file, gene_id_name_file):
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
        handle = Entrez.esearch(db="nuccore", term=gene_search)
        record = Entrez.read(handle)
        handle.close()
        id_list += record["IdList"][:1]  # set number of search results to read
    with open(tmp + gene_id_name_file, "w") as f:
        f.write("\n".join(id_list))

def ncbi_gi_to_genbank(id_list):
    # Get the genbank file of each gene by search for id list
    handle = Entrez.efetch(
        db="nucleotide",
        strand=1,  # plus if strand=1
        id=id_list,
        rettype="gbwithparts",
        retmode="text",
    )
    seq_record = handle.read()
    handle.close()
    return seq_record