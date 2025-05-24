# Code to Accelerate MSA computation using MMseqs2
# !mmseqs databases UniRef90 uniprotdb tmp
from Bio import SeqIO
from io import StringIO
import requests

def get_sequences(accession_numbers):
    base_url = "https://www.uniprot.org/uniprot/"
    sequences = {}

    for accession in accession_numbers:
        session = requests.Session()
        response = session.get(base_url + accession + ".fasta")
        session.close()

        sequence = next(SeqIO.parse(StringIO(response.text), "fasta"))
        sequences[accession] = str(sequence.seq)

    return sequences

import os, shutil
import pandas as pd  # assumed imported elsewhere

def find_homologs(df, column_name):
    results = []
    for index, row in df.iterrows():
        # 4.1 Create a per-row temp directory
        temp_dir = os.path.join("temp", str(index))
        os.makedirs(temp_dir, exist_ok=True)

        # 4.2 Write the query sequence to a FASTA file
        temp_fasta_file = os.path.join(temp_dir, f"temp_{index}.fasta")
        with open(temp_fasta_file, "w") as f:
            f.write(f">sequence_{index}\n{row[column_name]}\n")

        # 4.3 Define paths for MMseqs2 databases and results
        query_db  = os.path.join(temp_dir, "queryDB")
        result_db = os.path.join(temp_dir, "resultDB")
        result_m8 = os.path.join(temp_dir, "resultDB.m8")

        # 4.4 Run MMseqs2 commands
        os.system(f"mmseqs createdb {temp_fasta_file} {query_db} --dbtype 0")
        os.system(f"mmseqs search {query_db} uniprotdb {result_db} tmp --max-seqs 50000")
        os.system(f"mmseqs convertalis {query_db} uniprotdb {result_db} {result_m8}")

        # 4.5 Load and sort results, keep top 30 by bit-score
        cols = ['qId','tId','pctIdent','alnLen','mismatchCount',
                'gapOpenCount','qStart','qEnd','tStart','tEnd',
                'eVal','bitScore']
        result_df = pd.read_csv(result_m8, sep='\t', header=None, names=cols)
        accessions = ( result_df
                       .sort_values(by='bitScore', ascending=False)
                       ['tId']
                       .head(30)
                       .tolist() )

        # 4.6 Fetch the FASTA sequences for those top hits
        sequences = get_sequences(accessions)
        results.append(sequences)

        # 4.7 Clean up temp files
        shutil.rmtree(temp_dir)

    return results

if __name__ == "__main__":
    df = pd.read_csv("Data/full_type_seq.csv")  # Load your DataFrame
    all_homologs = find_homologs(df, "sequence")
    all_homologs_list = []
    for homolog in all_homologs:
        all_homologs_list.append(list(homolog.values()))
    df = df.assign(MSA=all_homologs_list)

    