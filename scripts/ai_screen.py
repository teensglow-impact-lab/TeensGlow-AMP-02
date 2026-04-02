"""
TeensGlow-AMP-02 | Module 3: AI-Driven Screening
=================================================
Screens bacteriocin candidates using ESM-2 protein language model embeddings
combined with a multi-dimensional TeensGlow Scoring System (TGS).

Scoring system (see Methods, Section 2.3):
    TGS = AMP_Potential × (TI/10 + SC/100) − Toxicity/3

where:
    AMP_Potential  : ESM-2 embedding-derived antimicrobial activity estimate
    TI             : Targeting Index for C. acnes (charge density + hydrophobicity)
    SC             : Sebum Compatibility Index (stability in lipid-rich sebaceous follicles)
    Toxicity       : Estimated hemolytic propensity (GRAVY-based + motif penalty)

Model: facebook/esm2_t6_8M_UR50D
    - 6 transformer layers, 20 attention heads, embedding dim = 320
    - Zero-shot inference (no task-specific fine-tuning)
    - Reference: Lin et al. (2023) Science 379:1123-1130

Output:
    results/ai_screening_results.csv
"""

import os
import torch
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from transformers import AutoTokenizer, EsmModel
from tqdm import tqdm


def compute_physicochemical(seq: str):
    """
    Compute key physicochemical descriptors for a peptide sequence.

    Returns
    -------
    charge : float
        Net charge at physiological pH (7.0).
    gravy : float
        Grand Average of hYdropathicity (GRAVY) index.
    mw : float
        Molecular weight (Da).
    ti : float
        Targeting Index — optimised for C. acnes penetration in sebaceous skin.
    sc : float
        Sebum Compatibility Index — stability in lipid-rich follicular environment.
    """
    analysis = ProteinAnalysis(seq)
    charge = analysis.charge_at_pH(7.0)
    gravy = analysis.gravy()
    mw = analysis.molecular_weight()

    # Targeting Index (TI): higher charge density + moderate hydrophobicity
    # favours electrostatic attraction to negatively charged C. acnes cell wall
    ti = (charge / len(seq)) * 8 + (gravy + 0.5) * 4

    # Sebum Compatibility (SC): reflects stability in lipid-rich microenvironment
    # Higher GRAVY → greater lipid-phase partitioning → better follicle penetration
    sc = float(np.clip((gravy + 1) * 35, 20, 100))

    return (
        round(charge, 2),
        round(gravy, 2),
        round(mw, 2),
        round(ti, 2),
        round(sc, 2)
    )


def screen_candidates(fasta_path: str, out_csv: str):
    """
    Score all candidate peptides using ESM-2 embeddings and the TGS formula.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file of candidate peptides.
    out_csv : str
        Output path for scored results CSV.

    Returns
    -------
    pd.DataFrame
        Sorted results table, highest GlowScore first.
    """
    model_name = "facebook/esm2_t6_8M_UR50D"
    print(f"[INFO] Loading ESM-2 model: {model_name}  (zero-shot inference mode)")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)
    model.eval()

    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"[INFO] Scoring {len(records)} candidate peptides...")

    results = []
    for rec in tqdm(records, desc="ESM-2 scoring"):
        seq = str(rec.seq)

        # --- Physicochemical descriptors ---
        charge, gravy, mw, ti, sc = compute_physicochemical(seq)

        # --- ESM-2 embedding (mean-pooled over residue dimension) ---
        inputs = tokenizer(seq, return_tensors="pt", add_special_tokens=False)
        with torch.no_grad():
            outputs = model(**inputs)
        embedding = outputs.last_hidden_state.mean(dim=1).squeeze().numpy()

        # --- AMP Potential: non-linear combination of embedding mean and charge ---
        amp_potential = np.tanh(np.mean(embedding) * 10 + (charge / 20)) * 50 + 50

        # --- Toxicity estimate: penalise high hydrophobicity and Leu-rich motifs ---
        # (Adolescent skin is typically more sensitive; lower irritancy is preferred)
        toxicity = (gravy + 0.5) * 15 + (20 if "LLL" in seq else 0)

        # --- TeensGlow GlowScore (TGS) ---
        glow_score = (amp_potential * (ti / 10 + sc / 100)) - (toxicity / 3)

        results.append({
            "ID": rec.id,
            "Sequence": seq,
            "Charge": charge,
            "GRAVY": gravy,
            "MW_Da": mw,
            "TargetingIndex": ti,
            "Sebum_Compat": sc,
            "AMP_Potential": round(amp_potential, 2),
            "Toxicity": round(toxicity, 2),
            "GlowScore": round(glow_score, 2),
            "Source": rec.description.split("Org:")[-1].strip()
                      if "Org:" in rec.description else "Unknown"
        })

    df = pd.DataFrame(results).sort_values(by="GlowScore", ascending=False)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)

    print(f"\n[DONE] Screening complete. Results saved to: {out_csv}")
    print(f"\n  Top 5 Candidates:")
    print(df[["ID", "GlowScore", "TargetingIndex", "Sebum_Compat", "Toxicity"]].head(5).to_string(index=False))
    return df


if __name__ == "__main__":
    screen_candidates(
        fasta_path="data/processed/candidates.fasta",
        out_csv="results/ai_screening_results.csv"
    )
