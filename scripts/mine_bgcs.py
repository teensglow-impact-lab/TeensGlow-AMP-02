"""
TeensGlow-AMP-02 | Module 2: BGC/ORF Mining
=============================================
Extracts small, cationically-charged open reading frames (ORFs) from
probiotic genomes as candidate bacteriocin sequences.

Filtering criteria (based on established AMP literature):
  - Length: 30-80 amino acids
  - Net charge at pH 7.0: >= +3.0
  - GRAVY index: -0.8 to +0.8
  - No low-complexity sequences (single AA > 50% composition)
  - Start codon: ATG only (strict)

Output:
  - data/processed/candidates.fasta  (FASTA, top 500 by charge)
  - data/processed/candidates.csv    (metadata table)
"""

import os
import gzip
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from tqdm import tqdm


def is_low_complexity(seq: str, threshold: float = 0.50) -> bool:
    """Return True if any single amino acid exceeds the composition threshold."""
    for aa in set(seq):
        if seq.count(aa) / len(seq) > threshold:
            return True
    return False


def find_candidate_orfs(nucleotide_seq, min_aa: int = 30, max_aa: int = 80):
    """
    Scan all 6 reading frames for ORFs meeting bacteriocin candidate criteria.

    Parameters
    ----------
    nucleotide_seq : Bio.Seq
        Full nucleotide sequence of a genomic contig.
    min_aa : int
        Minimum peptide length in amino acids.
    max_aa : int
        Maximum peptide length in amino acids.

    Returns
    -------
    list of dict
        Each dict contains peptide sequence and computed physicochemical properties.
    """
    candidates = []
    seq_str = str(nucleotide_seq).upper()
    reverse_complement = str(nucleotide_seq.reverse_complement()).upper()

    for strand, nuc in [(+1, seq_str), (-1, reverse_complement)]:
        for frame in range(3):
            # Use regex to find ATG...stop codon patterns efficiently
            for match in re.finditer(r"ATG(?:...)+?(?:TAA|TAG|TGA)", nuc[frame:]):
                orf_nuc = match.group()
                if min_aa * 3 <= len(orf_nuc) <= max_aa * 3:
                    protein = str(Seq(orf_nuc).translate(to_stop=True))
                    if len(protein) < min_aa:
                        continue
                    if is_low_complexity(protein):
                        continue

                    try:
                        analysis = ProteinAnalysis(protein)
                        charge = analysis.charge_at_pH(7.0)
                        gravy = analysis.gravy()

                        # Retain only cationic, moderately hydrophobic candidates
                        if charge >= 3.0 and -0.8 <= gravy <= 0.8:
                            candidates.append({
                                "Strand": strand,
                                "Peptide": protein,
                                "NetCharge": round(charge, 2),
                                "GRAVY": round(gravy, 2),
                                "Length": len(protein)
                            })
                    except Exception:
                        pass
    return candidates


def process_genomes(metadata_csv: str, out_fasta: str, max_candidates: int = 500):
    """
    Main pipeline: iterate over downloaded genomes, mine ORFs, filter,
    deduplicate, and export top candidates.
    """
    df = pd.read_csv(metadata_csv)
    all_candidates = []

    print(f"[INFO] Mining candidate ORFs from {len(df)} genome(s)...")
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Genomes processed"):
        local_path = row["LocalPath"]
        if not os.path.exists(local_path):
            print(f"  [WARNING] File not found, skipping: {local_path}")
            continue

        handle = gzip.open(local_path, "rt") if local_path.endswith(".gz") \
            else open(local_path, "r")

        for record in SeqIO.parse(handle, "fasta"):
            orfs = find_candidate_orfs(record.seq)
            for orf in orfs:
                orf["Organism"] = row["Organism"]
                orf["AssemblyID"] = row["ID"]
                all_candidates.append(orf)
        handle.close()

    candidates_df = pd.DataFrame(all_candidates)

    if candidates_df.empty:
        print("[WARNING] No candidates found. Consider relaxing filter thresholds.")
        return

    # Deduplicate by exact peptide sequence
    candidates_df = candidates_df.drop_duplicates(subset=["Peptide"])
    print(f"[INFO] Total unique candidates after deduplication: {len(candidates_df)}")

    # Rank by net charge (proxy for membrane interaction potential)
    if len(candidates_df) > max_candidates:
        print(f"[INFO] Retaining top {max_candidates} candidates by net charge.")
        candidates_df = candidates_df.sort_values(
            by="NetCharge", ascending=False
        ).head(max_candidates)

    # Export FASTA
    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
    records = [
        SeqRecord(
            Seq(row["Peptide"]),
            id=f"AMP_{i}_{row['AssemblyID']}",
            description=(
                f"Charge:{row['NetCharge']} "
                f"Gravy:{row['GRAVY']} "
                f"Org:{row['Organism']}"
            )
        )
        for i, (_, row) in enumerate(candidates_df.iterrows())
    ]
    SeqIO.write(records, out_fasta, "fasta")

    # Export metadata CSV
    csv_out = out_fasta.replace(".fasta", ".csv")
    candidates_df.to_csv(csv_out, index=False)
    print(f"[DONE] {len(candidates_df)} candidates saved to: {out_fasta}")


if __name__ == "__main__":
    process_genomes(
        metadata_csv="data/raw/assembly_metadata.csv",
        out_fasta="data/processed/candidates.fasta"
    )
