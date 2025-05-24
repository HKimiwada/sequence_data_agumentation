# One Sequence = 2 minutes 
import os
import time
import tempfile
import subprocess
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML

# === Settings ===
CSV_PATH   = "Data/full_type_seq.csv"
EMAIL      = "hiroki.kimiwada@keio.jp"   # NCBI requires a valid email
EVALUE     = 1e-3
TOP_K      = 100
OUTPUT_PKL = "data_with_msa.pkl"
MAX_WORKERS = 4                        # adjust to the number of CPU cores you can spare

Entrez.email = EMAIL

def get_top_hits(seq, evalue_thresh=EVALUE, top_k=TOP_K):
    result_handle = NCBIWWW.qblast("blastp", "nr", seq)
    record = NCBIXML.read(result_handle)
    hits = []
    for aln in record.alignments:
        if len(hits) >= top_k: break
        for hsp in aln.hsps:
            if hsp.expect < evalue_thresh:
                hits.append(aln.accession)
                break
    return hits

def fetch_seqs(accessions):
    records = []
    for acc in accessions:
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta")
        records.append(SeqIO.read(handle, "fasta"))
        handle.close()
    return records

def build_msa_for_sequence(idx, query_seq):
    """
    Full pipeline for one query:
      - BLAST → top hits
      - fetch hits → SeqRecords
      - add query sequence
      - MAFFT
      - parse alignment
    Returns: (idx, MultipleSeqAlignment)
    """
    # 1) BLAST
    hits = get_top_hits(query_seq)
    # 2) fetch and append query
    seqs = fetch_seqs(hits)
    seqs.append(SeqRecord(Seq(query_seq), id=f"query_{idx}"))
    # 3) write temp FASTA
    with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".fasta") as tmp_in:
        SeqIO.write(seqs, tmp_in.name, "fasta")
    tmp_out = tmp_in.name.replace(".fasta", "_aln.fasta")
    # 4) run MAFFT quietly
    subprocess.run(
        ["mafft", "--quiet", tmp_in.name],
        stdout=open(tmp_out, "w"),
        stderr=subprocess.DEVNULL
    )
    # 5) read alignment back
    alignment = AlignIO.read(tmp_out, "fasta")
    # cleanup
    os.remove(tmp_in.name)
    os.remove(tmp_out)
    return idx, alignment

def main():
    # load your data
    data = pd.read_csv(CSV_PATH)[["toughness", "sequence"]].head(2).dropna().reset_index(drop=True)
    n = len(data)
    msas = [None] * n

    # prepare the worker function
    worker = partial(build_msa_for_sequence)

    print(f"Spawning up to {MAX_WORKERS} workers to process {n} sequences…")
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(worker, idx, row["sequence"]): idx
            for idx, row in data.iterrows()
        }
        for i, future in enumerate(as_completed(futures), 1):
            idx, alignment = future.result()
            msas[idx] = alignment
            print(f"[{i:>3}/{n}] Done query_{idx}")
            # throttle to be kind to NCBI if you like:
            time.sleep(0.5)

    # attach and persist
    data["MSA"] = msas
    data.to_pickle(OUTPUT_PKL)
    print(f"\n✅ Finished. Pickled DataFrame with MSAs to {OUTPUT_PKL}")

if __name__ == "__main__":
    main()
