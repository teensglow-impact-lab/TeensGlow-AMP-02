"""
TeensGlow-AMP-02 | Module 4: Lead Peptide Truncation Optimisation
=================================================================
Performs systematic sliding-window truncation of the full-length lead
peptide to identify the minimal active core fragment.

Rationale: Shorter peptides (15-40 AA) typically:
  - Have lower synthesis cost (critical for commercialisation)
  - Show improved membrane selectivity
  - Are less susceptible to proteolytic degradation in sebum

Method:
  For each window size w (15 ≤ w ≤ 40):
    For each start position s:
      Extract subsequence, re-score with ESM-2 TGS formula

Output:
  results/alpha_optimization.csv   (all fragments + scores, sorted by GlowScore)
"""

import os
import torch
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from transformers import AutoTokenizer, EsmModel
from tqdm import tqdm


def compute_glow_score(seq: str, model, tokenizer) -> tuple:
    """
    Compute TeensGlow GlowScore (TGS) for a single peptide fragment.

    Returns
    -------
    glow_score : float
    targeting_index : float
    amp_potential : float
    """
    analysis = ProteinAnalysis(seq)
    charge = analysis.charge_at_pH(7.0)
    gravy = analysis.gravy()

    ti = (charge / len(seq)) * 8 + (gravy + 0.5) * 4
    sc = float(np.clip((gravy + 1) * 35, 20, 100))

    # ESM-2 embedding
    inputs = tokenizer(seq, return_tensors="pt", add_special_tokens=False)
    with torch.no_grad():
        outputs = model(**inputs)
    embedding = outputs.last_hidden_state.mean(dim=1).squeeze().numpy()

    # Length-adjusted AMP potential (peak activity ~15-30 AA for classical AMPs)
    len_weight = np.exp(-((len(seq) - 22) ** 2) / 300)
    amp_potential = np.tanh(np.mean(embedding) * 10 + (charge / 20)) * 50 + 50
    amp_potential *= len_weight

    toxicity = (gravy + 0.5) * 15
    glow_score = (amp_potential * (ti / 10 + sc / 100)) - (toxicity / 3)

    return round(glow_score, 2), round(ti, 2), round(amp_potential, 2)


def truncation_optimisation(
    lead_sequence: str,
    lead_id: str,
    out_csv: str,
    min_window: int = 15,
    max_window: int = 40
):
    """
    Exhaustive sliding-window scan to find the minimal active core.

    Parameters
    ----------
    lead_sequence : str
        Full amino acid sequence of the lead candidate.
    lead_id : str
        Identifier label for reporting purposes.
    out_csv : str
        Output path for the truncation results table.
    min_window : int
        Minimum fragment length to evaluate.
    max_window : int
        Maximum fragment length to evaluate.
    """
    print(f"[INFO] Loading ESM-2 for truncation optimisation of: {lead_id}")
    model_name = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)
    model.eval()

    results = []
    total = sum(
        len(lead_sequence) - w + 1
        for w in range(min_window, max_window + 1)
    )

    print(f"[INFO] Scanning {total} fragment(s) "
          f"(window: {min_window}–{max_window} AA)...")

    with tqdm(total=total, desc="Truncation scan") as pbar:
        for w in range(min_window, max_window + 1):
            for start in range(len(lead_sequence) - w + 1):
                fragment = lead_sequence[start:start + w]
                score, ti, pot = compute_glow_score(fragment, model, tokenizer)
                results.append({
                    "Fragment": fragment,
                    "Start_Pos": start + 1,
                    "End_Pos": start + w,
                    "Length_AA": w,
                    "GlowScore": score,
                    "TargetingIndex": ti,
                    "AMP_Potential": pot
                })
                pbar.update(1)

    df = pd.DataFrame(results).sort_values(by="GlowScore", ascending=False)
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    df.to_csv(out_csv, index=False)

    best = df.iloc[0]
    print(f"\n[DONE] Optimisation complete.")
    print(f"  Best fragment : {best['Fragment']}")
    print(f"  Position      : {int(best['Start_Pos'])}–{int(best['End_Pos'])}")
    print(f"  Length        : {int(best['Length_AA'])} AA")
    print(f"  GlowScore     : {best['GlowScore']}")
    print(f"  Results saved to: {out_csv}")


if __name__ == "__main__":
    # Alpha: top-ranked full-length lead from Module 3
    ALPHA_SEQUENCE = (
        "MRVALRLILKTQKKPCKYYMMRLLKLVIHLVKKLRLLLTRRHR"
        "IFIIERPTIMTLRAKRTRRQPLVIIICSY"
    )
    truncation_optimisation(
        lead_sequence=ALPHA_SEQUENCE,
        lead_id="Alpha",
        out_csv="results/alpha_optimization.csv"
    )
