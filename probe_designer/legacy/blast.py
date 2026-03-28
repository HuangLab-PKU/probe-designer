from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

def blast_online(tmp, tmp_output_pd, total_pre_binding_file_name, blast_results_file='5_blast_results.xml'):
    with open(tmp_output_pd, +total_pre_binding_file_name, 'r') as f:
        fasta_string = f.read()
    txid = [2697049]  # organism
    # Submit BLAST search and get handle object
    handle = NCBIWWW.qblast(program='blastn', megablast="yes",
                            database='refseq_rna', sequence=fasta_string,
                            url_base='https://blast.ncbi.nlm.nih.gov/Blast.cgi', 
                            format_object='Alignment',
                            format_type='Xml')
    with open(tmp+blast_results_file, 'w') as f:
        f.write(handle.read())


def blast_local(tmp, BDS_file_out_dir, total_pre_binding_file_name, blast_results_file='5_blast_results.xml', 
                 db="refseq_rna", task="megablast"): # db='refseq_rna' or 'refseq_select_rna'
    # Define the blastn command line
    blastn_cline = NcbiblastnCommandline(query=BDS_file_out_dir+total_pre_binding_file_name, db=db, task=task, 
                                        outfmt=5, out=tmp+blast_results_file) # outfmt=5 generates an XML output that includes alignment details
    
    # Execute the blast command
    stdout, stderr = blastn_cline()


def extract_blast_inf(tmp, blast_results_file, tmp_output_pd):
    # Extract interested information from blast_results
    align_num = []

    # read the id/plus-minus part/align_num
    with open(tmp+blast_results_file, 'r') as blast_output:
        blast_records = NCBIXML.parse(blast_output)
        loca = 0
        for blast_record in blast_records:
            align_accession = []
            align_descrip_list = []
            # get align num of each binding site
            length = len(blast_record.alignments)
            align_num.append(length)
            for i in range(length):
                descrip = blast_record.descriptions[i].title.split('|')
                # get accession and descrip of each align seq
                align_accession.append(descrip[3])
                align_descrip_list.append(descrip[-1])
            tmp_output_pd.loc[loca, 'align_accession'] = '|'.join(str(_) for _ in align_accession)
            
            # add align_descrip to df
            tmp_output_pd.loc[loca, 'align_descrip'] = '|'.join(str(_) for _ in align_descrip_list)
            
            # get plus/minus of each align seq
            p_m = [blast_record.alignments[_].hsps[0].frame[1] for _ in range(length)]
            
            # add plus/minus to df
            try: tmp_output_pd.loc[loca, 'plus/minus'] = ','.join([str(_) for _ in p_m])
            except: tmp_output_pd.loc[loca, 'plus/minus'] = 'NAN'
            
            loca += 1

    tmp_output_pd['align_num'][:len(align_num)] = align_num
    
    return tmp_output_pd