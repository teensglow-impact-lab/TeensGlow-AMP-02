import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_benchmarks(data_csv, out_path):
    df = pd.read_csv(data_csv)
    
    # Separate data: Candidates vs Benchmarks
    df["Category"] = df["ID"].apply(lambda x: "Benchmark" if "Benchmark" in x else "TeensGlow Candidate")
    
    plt.figure(figsize=(10, 7), dpi=300)
    sns.set_style("darkgrid")
    
    # Plot scatter: Targeting Index vs AMP Potential
    scatter = sns.scatterplot(
        data=df, 
        x="TargetingIndex", 
        y="AMP_Potential", 
        hue="Category", 
        size="Charge", 
        sizes=(20, 200),
        palette={"Benchmark": "red", "TeensGlow Candidate": "dodgerblue"},
        alpha=0.7
    )
    
    # Annotate Top-1 (Alpha-Core) and reference Benchmarks
    top_1 = df.iloc[0]
    plt.text(top_1["TargetingIndex"]+0.1, top_1["AMP_Potential"], "Alpha", fontsize=10, weight='bold')
    
    for i, row in df[df["Category"] == "Benchmark"].iterrows():
        plt.text(row["TargetingIndex"]+0.1, row["AMP_Potential"], row["ID"].split("_")[-1], fontsize=9)

    plt.title("Comparative Analysis: TeensGlow Leads vs Reference AMPs", size=15)
    plt.xlabel("Targeting Index (C. acnes specificity)", size=12)
    plt.ylabel("Antimicrobial Potential (AI embedding score)", size=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    
    plt.savefig(out_path, bbox_inches="tight")
    print(f"[+] Comparative scatter plot saved to {out_path}")

if __name__ == "__main__":
    if not os.path.exists("results/figures"): os.makedirs("results/figures")
    plot_benchmarks("results/polished_sci_data.csv", "results/figures/comparative_plot.png")
