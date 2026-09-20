import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")

ADJ_P_THR = 0.05
LFC_THR = 0.585

gene_sets = {}
for org, fname in [("Kidney", "kidney_DEGs.csv"),
                  ("Liver",  "liver_DEGs.csv"),
                  ("Lung",   "lung_DEGs.csv"),
                  ("Skin",   "skin_DEGs.csv")]:
    df = pd.read_csv(os.path.join(RESULTS, fname))
    sig = df[(df["adj_p_value"] < ADJ_P_THR) & (df["logFC"].abs() >= LFC_THR)].copy()
    sig["gene"] = sig["gene"].astype(str).str.strip().str.upper()
    gene_sets[org] = set(sig["gene"].unique())

print("Loaded organ sets:")
for org, s in gene_sets.items():
    print(f"  {org}: {len(s):,} genes")

core = gene_sets["Kidney"] & gene_sets["Liver"] & gene_sets["Lung"] & gene_sets["Skin"]
print(f"4-organ intersection: {len(core)} genes")

upset_data = from_contents(gene_sets)

fig = plt.figure(figsize=(14, 8), dpi=300)
upset = UpSet(upset_data, subset_size="count", show_counts=True,
              sort_by="cardinality", show_percentages=True,
              element_size=42)
upset.plot(fig=fig)

fig.suptitle(
    f"Overlap of Significant Fibrosis DEGs Across 4 Organs (Corrected Discovery: N=14 Cohorts, 1,069 Samples)\n"
    f"[Kidney: {len(gene_sets['Kidney']):,} | Liver: {len(gene_sets['Liver']):,} | Lung: {len(gene_sets['Lung']):,} | Skin: {len(gene_sets['Skin']):,} | 4-Organ Core: {len(core)}]",
    fontsize=13, fontweight="bold", y=1.03, color="#1a1a4a"
)

out_path = os.path.join(BASE, "plots", "upset_plot_4organs_corrected.png")
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"[SAVED] {out_path}")
plt.close(fig)
