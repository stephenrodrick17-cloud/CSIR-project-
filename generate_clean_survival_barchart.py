import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(BASE, "results", "final_master_evidence_table.csv")

df = pd.read_csv(CSV_PATH)

# Parse fractions like "2/4" -> 2
def parse_fraction(val):
    if pd.isna(val):
        return 0
    s = str(val).strip()
    if "/" in s:
        return int(s.split("/")[0])
    return int(s)

df["val1_sig_num"] = df["val1_significant_organs"].apply(parse_fraction)
df["val2_sig_num"] = df["val2_significant_organs"].apply(parse_fraction)
df["val2_conc_num"] = df["val2_concordant_organs"].apply(parse_fraction)

# Sort by tier, then val2_sig_num desc, then val1_sig_num desc, then ml_diagnostic_auc desc
tier_order = {"Tier 1 (Full Spectrum)": 1, "Tier 2 (Validation-Only)": 2, "Not Supported": 3}
df["tier_rank"] = df["final_evidence_tier"].map(tier_order)
df = df.sort_values(by=["tier_rank", "val2_sig_num", "val1_sig_num", "ml_diagnostic_auc"], ascending=[True, False, False, False]).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(16, 8), dpi=300)

x = np.arange(len(df))
width = 0.38

rects1 = ax.bar(x - width/2, df["val1_sig_num"], width, label="Val 1 Significant Organs (FDR < 0.05, /4)", color="#2b5c8f", alpha=0.9, edgecolor="black", linewidth=0.8)
rects2 = ax.bar(x + width/2, df["val2_sig_num"], width, label="Val 2 Held-Out Significant Organs (FDR < 0.05, /4)", color="#d95f02", alpha=0.9, edgecolor="black", linewidth=0.8)

ax.set_ylabel("Number of Statistically Significant Organs (out of 4)", fontsize=13, fontweight="bold")
ax.set_title("Validation Survival Analysis of 24 Clean Discovery Core ECM Genes\nStratified by Validation 1 (Internal Cohorts) vs. Validation 2 (Held-Out Cohorts)", fontsize=15, fontweight="bold", pad=15)
ax.set_xticks(x)
ax.set_xticklabels(df["gene"], rotation=45, ha="right", fontsize=11, fontweight="bold")
ax.set_ylim(0, 5.2)
ax.axhline(2, color="gray", linestyle="--", linewidth=1.2, alpha=0.7, label="Multi-Organ Replicability Threshold (>= 2/4)")
ax.axhline(4, color="forestgreen", linestyle=":", linewidth=1.5, alpha=0.8, label="Universal 4/4 Organ Replication")

# Color code x-tick labels by tier
for i, tick in enumerate(ax.get_xticklabels()):
    tier = df.loc[i, "final_evidence_tier"]
    if tier == "Tier 1 (Full Spectrum)":
        tick.set_color("#006400") # dark green
    elif tier == "Tier 2 (Validation-Only)":
        tick.set_color("#1a365d") # dark navy
    else:
        tick.set_color("#8b0000") # dark red

# Add value labels on top of bars
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        if height > 0:
            ax.annotate(f"{int(height)}",
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha="center", va="bottom", fontsize=10, fontweight="bold")

autolabel(rects1)
autolabel(rects2)

# Annotate tiers
tier1_end = len(df[df["tier_rank"] == 1]) - 0.5
tier2_end = tier1_end + len(df[df["tier_rank"] == 2])

ax.axvline(tier1_end, color="#006400", linestyle="-.", linewidth=1.5, alpha=0.7)
ax.axvline(tier2_end, color="#8b0000", linestyle="-.", linewidth=1.5, alpha=0.7)

ax.text(tier1_end / 2, 4.75, "TIER 1: Full Spectrum (N=4)\n(Val1 + Val2 >=2/4 + ML Stability >=4/5)", 
        ha="center", va="center", fontsize=9, fontweight="bold", color="#004d00",
        bbox=dict(boxstyle="round,pad=0.3", fc="#e6f4ea", ec="#006400", lw=1))

ax.text(tier1_end + (tier2_end - tier1_end) / 2, 4.75, "TIER 2: Biological Validation-Only (N=14)\n(Val1 4/4 + Val2 >=2/4 Replicated)", 
        ha="center", va="center", fontsize=9, fontweight="bold", color="#1a365d",
        bbox=dict(boxstyle="round,pad=0.3", fc="#e8f0fe", ec="#1a365d", lw=1))

ax.text(tier2_end + (len(df) - tier2_end) / 2, 4.75, "NOT SUPPORTED (N=6)\n(Val2 Replicated in 1/4 organ)", 
        ha="center", va="center", fontsize=9, fontweight="bold", color="#8b0000",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fce8e6", ec="#8b0000", lw=1))

ax.grid(axis="y", linestyle=":", alpha=0.6)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=10.5, framealpha=0.95)

plt.tight_layout()

out_path = os.path.join(BASE, "plots", "clean_ecm_validation_survival_barchart.png")
plt.savefig(out_path, dpi=300)
print(f"[SAVED] {out_path}")
plt.close()
