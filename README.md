# probe_designer

This program is intended to design padlock probes for **spatial transcriptome** or **multiplex RNA FISH**. When input with gene list you are interested in, the program will search for the corresponding mRNA sequence for the gene and, with some strategy, search for binding sites of specific length you assign.

## Quickstart
*Because the package for RNA 2nd structure detection is linux only, this code is recommanded to run on **linux**, **mac os** or **windows subsystem for linux**.*

### Step 1: create the environment

Clone the repo to your app/repo dir:
```bash
git clone --depth 1 https://github.com/tangmc0210/probe_designer
```

Create the environment using conda:
```bash
conda env create --file environment.yml
```

### Step 2: prepare your gene list
Create a directory for you task, like `exmple_dataset`. Prepare a .xlsx file including column "gene_name" which include target genes as `data/gene_list.xlsx`.

| gene_name | ... |
| --------- | --- |
| CD3D      | ... |
| CD4       | ... |
| CD8A      | ... |
| ...       | ... |


### Step 3: select binding sites 
Run `binding_search_ensmbl.ipynb` step by step. You need to adjust the workpath, the organism type and you may also revise the gene_list loading part based on your input file content. Note that because the NCBIWWW Blast is not stable when uploading big files, you may need to blast the fasta on [NCBI BLAST webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi) using the file generated in your result path: `results/datetime/total_bds_candidate_blast.fasta`, download the results in `XML` format and put it as `results/datetime/blast_results`.

### Step 4: check not found genes

In this step, the program identifies genes from the provided list that were not found in the binding site search results. It reads the gene names from `marker_gene_list.xlsx`, merges the results from `probes_wanted.xlsx` files, and standardizes the gene names for comparison. The program then filters for any genes that do not appear in the results and records them in `to_search.txt` for further review. This ensures that any unanalyzed genes are documented for follow-up.

### Step 5: iterate until getting all marker genes
Using binding_searcher files of ensembl or NCBI to get to search genes, you can loosen the conditions or manual select proper genes in ncbi web database. As an example, you can go to [ncbi webpage](https://www.ncbi.nlm.nih.gov/) and search in 'nucleotide' database like: `Apod, mouse, mRNA`. And record the corresponding entry's GI in `binding_searcher_NCBI.ipynb` file tree.


> ## Quick start with `main.py` 
> ***To develop***
> 
> ``` 
> workdir = '~/probe_designer/dataset'
> ```
> 
> ### Run with interactive command window
>
> Run code in proper environment:
> 
> ```powershell
> python main.py
> ```
> 
> And prepare files following commands.

## Binding site search strategy

The binding site search strategy has been improved by optimizing the search for binding positions that are maximally distant from each other. The program scans every possible point along the mRNA sequence and selects binding positions based on specific criteria, enhancing the accuracy of the search process.

The selection of binding sites follows two sets of rules: pre-BLAST rules, which filter sequences before conducting the BLAST search, and post-BLAST rules, which refine the selection based on sequence specificity after BLAST analysis. Pre-BLAST rules include criteria such as G content and melting temperature, while post-BLAST rules ensure the binding sites are specific to the target gene in the relevant organisms.

## Detailed description of ensenbl pipeline
### File Tree
```bash
RUN_ID
└──results
   └──20240910_134726_ensembl
     │   sequence_of_all.json
     │   shortest_isoforms.json
     │   bds_candidate_num.json
     │   total_bds_candidate.fasta
     │   total_bds_candidate_blast.fasta
     │   blast_results.xml    # manual
     │   probes_candidates.xlsx
     │   probes_wanted.xlsx
     │
     └───bds_candidate
            *.fasta
```
### Get Gene Information from Ensembl Dataset

The program starts by retrieving gene information from the Ensembl database based on a provided list of gene names. It utilizes the Ensembl API to gather relevant sequences, ensuring the sequences correspond to the specified organism. Note that gene names may require formatting adjustments depending on the species (e.g., uppercase for human genes).

This step generates the file `results/datetime/sequence_of_all.json`, which contains all retrieved sequences for further analysis.

### Download mRNA Sequences

Using the gene information obtained, the program fetches the mRNA sequences from Ensembl. If any retrieval attempts fail, the program retries up to three times to ensure robustness. This approach allows for the collection of the shortest isoform for each gene, optimizing the binding site search process.

This step produces the file `results/datetime/shortest_isoforms.json`, which lists the shortest isoforms for all genes of interest.

### Binding Site Searcher

As outlined in the **Binding Site Search Strategy**, the binding site search is optimized by locating binding positions that are maximally distant from each other. The program scans the mRNA sequences and identifies potential binding sites based on specific parameters such as binding site length and G content.

This section significantly improves hybridization performance and results in multiple FASTA files saved in `results/datetime/bds_candidate/*.fasta`. The main output file, `total_bds_candidate.fasta`, is prepared for subsequent BLAST analysis, and a summary of valid sequences is saved in `results/datetime/bds_candidate_num.json`.

### BLAST and Decoding of Results

After obtaining the candidate binding sites, you can perform a BLAST search on `total_bds_candidate_blast.fasta` to produce `results/datetime/blast_results.xml`. This file contains the alignment results for further interpretation. You may choose to run a local BLAST or utilize the NCBI BLAST web service for this task.

1. **Tip:** The web BLAST has a limit of about 1200 sequences per submission; consider splitting your input if you have many target genes.
2. **Local BLAST Setup:** For local BLAST, refer to [Download BLAST Software and Databases](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html).

### Select Wanted Binding Sites by Specific Rules

The program employs a two-tiered selection process to identify desired binding sites: pre-BLAST rules and post-BLAST rules. Pre-BLAST rules are applied before the BLAST analysis to filter out unsuitable sequences, while post-BLAST rules refine the selection based on alignment results.

#### Current Pre-BLAST Rules:
1. G content is between 40% and 70%.
2. Consecutive Gs do not exceed 4.
3. Melting temperature (Tm) of both arms of binding sites is between 45–65 °C.
4. Maximum free energy of RNA secondary structure is above -10 kcal/mol.
5. The sequence type must be mRNA and coding sequence (as determined by the GenBank file).

#### Current Post-BLAST Rules:
1. Binding site sequences must be specific to the target gene in relevant organisms (determined by BLAST results).
2. Binding site sequences should be positioned upstream or downstream of the target gene.

### Merge the Results

Run `merge.ipynb` to consolidate the results obtained from the above steps. This notebook compiles the data into an Excel file detailing the selected binding sites for each gene, along with a record of any missing genes in the pipeline, saved as `results/_tosearch.txt`.


## Detailed description of NCBI pipeline
### File Tree
```bash
RUN_ID
└──results
   └───20240910_140515_NCBI
      │   gene_name_list.txt    # manual
      │   gene_id_list.txt      # manual
      │   gene_seq_in_file.gb
      │   bds_candidate_num.json
      │   total_bds_candidate.fasta
      │   total_bds_candidate_blast.fasta
      │   blast_results.xml     # manual
      │   probes_candidates.xlsx
      │   probes_wanted.xlsx
      │
      └───bds_candidate
               SOX10_bds_candidate.fasta
```
### Get gene_id from ncbi dataset

Given the gene name list, the program search for target sequence using entrez provided by ncbi. However, the search results are usually mismatched so I recommand to search on website and choose target GI for next step analysis.

This step returns file `results/datetime/gene_id_list.txt`.

### Get genbank file of each mRNA from ncbi dataset

The program downloads rna seq from ncbi according to GI number given in tmp file "1_id_list.txt". So prepare the id_list file yourself by searching on web if you find it awful using auto search in previous step.

This step returns file `results/datetime/gene_seq_in_file.gb`

### Binding site searcher

As mentioned in **Binding site search strategy**, current search stratagy is bruteforce. The program will begin in the middle of a mRNA sequence and extend with a step length(length of mRNA sequence divided by binding site length) to both side.

This part is most likely to be improved if you want to get better hybriding performance.

This step returns a series of .fasta file as `results/datetime/bds_candidate/*.fasta` where `total_bds_candidates.fasta` is to blast. And it returns a `results/datetime/tmp/pre_binding_num.json` file which shows valid seq got from box searcher.

### Blast and decoding of blast_results

You can perform local blast or web blast on `total_bds_candidate_blast.fasta` to get `results/datetime/blast_results.xml` file for decoding. You will need to build your local blast if you need to perform blast in this code. Otherwise, you can get blast results by uploading the file to blast web and download the blast results as `results/datetime/blast_results.xml`.

1. Tip1: web blast length is limited to about 1200 sequences once, so split the gene name list if target gene number is too large.
2. How to build local blast: [Download BLAST Software and Databases](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

### Select wanted binding site by specific rules

Current rules are divided into two kinds: pre_blast rules and post_blast rules. Pre_blast rules work before blast the potential binding site sequences and can often exempt most of sequences while post_blast rules work after blast and are intended to get sequence specifity information.

Current before-blast rules:

1. G percentage is between 40% and 70%;
2. G is not consecutive up to 5 or more;
3. Tm of both arms of bds are between 45-65 °C;
4. Max free energy of RNA 2nd structure is above -10 kcal/mol; 
5. Sequence type is mRNA and coding sequence (judge by genbank file);

Current post-blast rules:

1. Binding site sequence is specific to this gene in interested organisms (judge by blast results);
2. Binding site sequence is plus/minus to target gene;

### Merge the results got from above

run `merge.ipynb` to merge results got above and return a xlsx file of wanted binding sites for each gene. And return the missing genes in the pipeline. `results/_tosearch.txt`.

## Perspective and plans

### Improve the searching strategy
- [x] RNA 2nd structure detection. 2024.6 done.
- [x] Most far away bds: The searching stategy is optimized by searching every possible points and select the binding positions most far away. 2023.6 done.

