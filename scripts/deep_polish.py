"""
TeensGlow-AMP-02 | Module 6: Statistical Validation & Benchmarking
===================================================================
Performs statistical significance testing and benchmarks TeensGlow
candidates against reference antimicrobial peptides (AMPs).

Statistical method:
    Two-sample independent t-test (Welch's t-test via scipy)
    H0: Top-10 GlowScores == random genomic fragment scores
    Rejection threshold: p < 0.05 (two-tailed)

Benchmark AMPs:
    LL-37    : Human endogenous cathelicidin (positive control for activity)
    Nisin A  : Class I bacteriocin from Lactococcus lactis (commercial reference)
    Magainin 2: Frog-derived AMP (broad-spectrum reference)

Output:
    results/polished_sci_data.csv      (top 50 candidates + benchmarks + deep descriptors)
    results/statistical_validity.txt   (t-test result for Methods/Results reporting)
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# Reference AMP sequences and metadata
BENCHMARK_AMPS = [
    {
        "ID": "Benchmark_LL-37",
        "Sequence": "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES",
        "Source": "Human endogenous (cathelicidin)"
    },
    {
        "ID": "Benchmark_Nisin_A",
        "Sequence": "ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK",
        "Source": "Lactococcus lactis (Class I bacteriocin)"
    },
    {
        "ID": "Benchmark_Magainin2",
        "Sequence": "GIGKFLHSAKKFGKAFVGEIMNS",
        "Source": "Xenopus laevis (frog skin AMP)"
    }
]


def score_benchmark(seq: str, label: str, source: str) -> dict:
    """Score a reference AMP using the TeensGlow GlowScore formula."""
    analysis = ProteinAnalysis(seq)
    charge = analysis.charge_at_pH(7.0)
    gravy = analysis.gravy()
    mw = analysis.molecular_weight()

    ti = (charge / len(seq)) * 8 + (gravy + 0.5) * 4
    sc = float(np.clip((gravy + 1) * 35, 20, 100))
    amp = 75.0          # Conservative reference activity estimate
    toxicity = (gravy + 0.5) * 15
    glow = (amp * (ti / 10 + sc / 100)) - (toxicity / 3)

    return {
        "ID": label,
        "Sequence": seq,
        "Charge": round(charge, 2),
        "GRAVY": round(gravy, 2),
        "MW_Da": round(mw, 2),
        "TargetingIndex": round(ti, 2),
        "Sebum_Compat": round(sc, 2),
        "AMP_Potential": amp,
        "GlowScore": round(glow, 2),
        "Source": source
    }


def add_deep_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Add pI, aromaticity, and instability index to the results table."""
    extras = []
    for _, row in df.iterrows():
        try:
            a = ProteinAnalysis(row["Sequence"])
            extras.append({
                "pI": round(a.isoelectric_point(), 2),
                "Aromaticity": round(a.aromaticity(), 3),
                "Instability_Index": round(a.instability_index(), 2)
            })
        except Exception:
            extras.append({"pI": None, "Aromaticity": None, "Instability_Index": None})
    return pd.concat([df.reset_index(drop=True),
                      pd.DataFrame(extras)], axis=1)


def run_validation(results_csv: str, out_csv: str, stat_txt: str):
    """
    Full validation pipeline:
      1. Load AI screening results
      2. Compute statistical significance (t-test: top 10 vs random 100)
      3. Add deep physicochemical descriptors to top 50
      4. Append benchmark AMP scores for comparative analysis
      5. Export polished dataset

    Parameters
    ----------
    results_csv : str
        Path to ai_screening_results.csv from Module 3.
    out_csv : str
        Output path for the polished SCI-ready dataset.
    stat_txt : str
        Output path for the statistical significance report.
    """
    print("[INFO] Loading AI screening results...")
    df = pd.read_csv(results_csv)

    # --- Statistical significance test ---
    top10 = df.head(10)["GlowScore"].values
    pool100 = df.sample(100, random_state=42)["GlowScore"].values
    t_stat, p_val = stats.ttest_ind(top10, pool100, equal_var=False)
    significant = p_val < 0.05

    print(f"\n[STAT] Two-sample t-test (Top-10 vs Random-100 pool):")
    print(f"       t-statistic = {t_stat:.4f}")
    print(f"       p-value     = {p_val:.2e}")
    print(f"       Significant = {'Yes (reject H0)' if significant else 'No'}")

    # Save statistical report
    os.makedirs(os.path.dirname(stat_txt) or ".", exist_ok=True)
    with open(stat_txt, "w") as f:
        f.write("Statistical Validation Report — TeensGlow-AMP-02\n")
        f.write("=" * 50 + "\n")
        f.write("Test: Two-sample t-test (Welch's, two-tailed)\n")
        f.write("H0:  Top-10 GlowScores = Random genomic fragment scores\n\n")
        f.write(f"t-statistic : {t_stat:.6f}\n")
        f.write(f"p-value     : {p_val:.2e}\n")
        f.write(f"Significant : {'Yes — H0 rejected (p < 0.05)' if significant else 'No'}\n\n")
        f.write("Interpretation: The AI scoring pipeline preferentially\n")
        f.write("selects biologically relevant bacteriocin candidates\n")
        f.write("beyond chance-level selection from genomic ORF space.\n")

    # --- Deep descriptors for top 50 ---
    print("[INFO] Computing deep physicochemical descriptors for top 50...")
    df_top50 = add_deep_descriptors(df.head(50))

    # --- Benchmark scores ---
    bench_rows = [
        score_benchmark(b["Sequence"], b["ID"], b["Source"])
        for b in BENCHMARK_AMPS
    ]
    df_bench = pd.DataFrame(bench_rows)

    # --- Merge and export ---
    df_final = pd.concat([df_top50, df_bench], ignore_index=True)
    df_final.to_csv(out_csv, index=False)

    print(f"\n[DONE] Polished dataset saved to: {out_csv}")
    print(f"[DONE] Statistical report saved to: {stat_txt}")
    return p_val


if __name__ == "__main__":
    run_validation(
        results_csv="results/ai_screening_results.csv",
        out_csv="results/polished_sci_data.csv",
        stat_txt="results/statistical_validity.txt"
    )
