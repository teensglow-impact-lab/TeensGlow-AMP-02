import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from math import pi

def make_radar_chart(df, out_path):
    # Select top 3 candidates for visualization
    top_3 = df.head(3)
    
    # Define radar chart dimensions (Activity, Targeting, Safety, Sebum Compatibility)
    categories = ["Activity", "Targeting", "Safety", "SebumCompat"]
    N = len(categories)
    
    # Data normalization to 0-100 scale for visual consistency
    def normalize(column):
        if column.max() == column.min(): return column * 0 + 50
        return (column - column.min()) / (column.max() - column.min()) * 60 + 40

    plot_data = pd.DataFrame({
        "Label": top_3["ID"],
        "Activity": normalize(top_3["AMP_Potential"]),
        "Targeting": normalize(top_3["TargetingIndex"]),
        "Safety": normalize(100 - top_3["Toxicity"]),
        "SebumCompat": normalize(top_3["Sebum_Compat"]) 
    })

    # Set up polar coordinate mapping for radar chart
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
    
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    
    # Background and grid styling properties
    plt.xticks(angles[:-1], categories, color="gray", size=12)
    ax.set_rlabel_position(0)
    plt.yticks([25, 50, 75, 100], ["25", "50", "75", "100"], color="grey", size=10)
    plt.ylim(0, 100)
    
    colors = ["#ff595e", "#1982c4", "#8ac926"]
    
    for i, (idx, row) in enumerate(plot_data.iterrows()):
        values = row[categories].values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2, linestyle="solid", label=row["Label"], color=colors[i])
        ax.fill(angles, values, color=colors[i], alpha=0.1)
        
    plt.legend(loc="upper right", bbox_to_anchor=(0.1, 0.1))
    plt.title("TeensGlow-AMP-02: Adolescent Acne Specificity Radar Chart", size=16, y=1.1)
    
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"[+] Radar chart saved to {out_path}")

if __name__ == "__main__":
    try:
        csv_path = "results/ai_screening_results.csv"
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            make_radar_chart(df, "results/figures/radar_chart.png")
        else:
            print("[!] Screening results CSV not found.")
    except Exception as e:
        print(f"[!] Plotting failed: {e}")
