# -*- coding: utf-8 -*-
"""
Module 6: GO / KEGG / Disease Ontology Enrichment Analysis
===========================================================
Tool   : gseapy v1.3.1 — gseapy.enrichr (live Enrichr REST API)
API    : https://maayanlab.cloud/Enrichr/
Stats  : Fisher's exact test with BH-FDR correction (performed server-side by Enrichr)

GENE SET 1 (9 genes) : COL15A1, COL1A1, COL3A1, SERPINE2, SERPINF2, LAMC3, LTBP2, MDK, SVEP1
GENE SET 2 (5 genes) : COL15A1, COL1A1, COL3A1, SERPINE2, SERPINF2

Libraries queried (same for both sets):
  1. GO_Biological_Process_2023
  2. GO_Molecular_Function_2023
  3. GO_Cellular_Component_2023
  4. KEGG_2021_Human
  5. DisGeNET

No genes are added, substituted, or bridged beyond the two locked sets above.
Any gene symbol that fails to map in Enrichr will be reported explicitly.

Outputs (SEPARATE — never merged):
  Gene Set 1 (9-gene):
    results/go_bp_9genes.csv
    results/go_mf_9genes.csv
    results/go_cc_9genes.csv
    results/kegg_9genes.csv
    results/do_9genes.csv
    plots/enrichment_9genes.png

  Gene Set 2 (5-gene):
    results/go_bp_5genes.csv
    results/go_mf_5genes.csv
    results/go_cc_5genes.csv
    results/kegg_5genes.csv
    results/do_5genes.csv
    plots/enrichment_5genes.png
"""

import os, sys, time, datetime

# Force UTF-8 stdout on Windows to avoid cp1252 encoding errors
if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout = open(sys.stdout.fileno(), mode="w", encoding="utf-8", buffering=1)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gseapy as gp

os.makedirs("results", exist_ok=True)
os.makedirs("plots",   exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# AUDIT HEADER
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("GO / KEGG / Disease Ontology Enrichment — LIVE Enrichr API")
print(f"gseapy version : {gp.__version__}")
print(f"Run timestamp  : {datetime.datetime.utcnow().isoformat()} UTC")
print("=" * 72)

# ─────────────────────────────────────────────────────────────────────────────
# LOCKED GENE SETS  (no additions, no substitutions)
# ─────────────────────────────────────────────────────────────────────────────
GENE_SETS = {
    "9genes": {
        "label": "9-Gene Consensus Hub Set",
        "genes": ["COL15A1", "COL1A1", "COL3A1",
                  "SERPINE2", "SERPINF2",
                  "LAMC3", "LTBP2", "MDK", "SVEP1"],
    },
    "5genes": {
        "label": "5-Gene Unanimous Hub Set",
        "genes": ["COL15A1", "COL1A1", "COL3A1",
                  "SERPINE2", "SERPINF2"],
    },
}

# Enrichr library → short tag → file suffix
LIBRARIES = [
    ("GO_Biological_Process_2023", "GO_BP",  "go_bp"),
    ("GO_Molecular_Function_2023", "GO_MF",  "go_mf"),
    ("GO_Cellular_Component_2023", "GO_CC",  "go_cc"),
    ("KEGG_2021_Human",            "KEGG",   "kegg"),
    ("DisGeNET",                   "DO",     "do"),
]

# Colour palette for plots (one colour per library)
LIB_COLORS = {
    "GO_BP": "#2563eb",
    "GO_MF": "#0891b2",
    "GO_CC": "#059669",
    "KEGG":  "#d97706",
    "DO":    "#dc2626",
}

TOP_N = 15          # terms per panel in the dot-plot figure
FDR_THRESH = 0.05   # significance line drawn on plots


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: run one Enrichr query, return clean DataFrame
# ─────────────────────────────────────────────────────────────────────────────
def run_enrichr(gene_list, library, set_tag, lib_tag):
    """
    Calls gseapy.enrichr (live REST API).
    Returns the full, unfiltered results DataFrame.
    """
    print(f"\n  → [LIVE API CALL] gseapy.enrichr("
          f"gene_list={gene_list}, "
          f"gene_sets='{library}', "
          f"organism='human', outdir=None, verbose=True)")
    print(f"     Library : {library}")
    print(f"     Gene set: {set_tag}  ({len(gene_list)} genes)")
    print(f"     Genes   : {', '.join(gene_list)}")

    t0 = time.time()
    enr = gp.enrichr(
        gene_list   = gene_list,
        gene_sets   = library,
        organism    = "human",
        outdir      = None,
        verbose     = True,
    )
    elapsed = time.time() - t0

    df = enr.results.copy()
    n_terms = len(df)
    print(f"     ✓ Returned {n_terms} terms in {elapsed:.1f}s")

    if df.empty:
        print(f"     ⚠  No terms returned for {library} / {set_tag}")
        return df

    # Rename columns to a consistent schema
    df = df.rename(columns={
        "Term"          : "Term",
        "Overlap"       : "Overlap",
        "P-value"       : "P_value",
        "Adjusted P-value": "Adj_P_value",
        "Old P-value"   : "Old_P_value",
        "Old Adjusted P-value": "Old_Adj_P_value",
        "Odds Ratio"    : "Odds_Ratio",
        "Combined Score": "Combined_Score",
        "Genes"         : "Genes",
    })

    # Add provenance columns
    df.insert(0, "Library",  library)
    df.insert(1, "Lib_Tag",  lib_tag)
    df.insert(2, "Gene_Set", set_tag)

    # Explicit: list any submitted genes absent from the returned "Genes" column
    all_returned_genes = set()
    for g_str in df["Genes"].dropna():
        all_returned_genes.update([g.strip().upper() for g in g_str.split(";")])

    missing = [g for g in gene_list if g.upper() not in all_returned_genes]
    if missing:
        print(f"     ⚠  Gene(s) NOT appearing in any returned term: {missing}")
        print(f"        (May indicate an alias mismatch in Enrichr's background)")
    else:
        print(f"     ✓ All {len(gene_list)} submitted genes appear in at least one returned term.")

    return df


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: dot-plot figure for one gene set (all 5 libraries, top N each)
# ─────────────────────────────────────────────────────────────────────────────
def make_dotplot(results_dict, set_tag, set_label, out_png, top_n=TOP_N):
    """
    Builds a 5-panel dot-plot figure (one panel per library).
    X-axis = -log10(Adj_P_value), dot size ∝ overlap count,
    dot colour = library colour.
    """
    n_panels = len(LIBRARIES)
    fig = plt.figure(figsize=(18, 5 * n_panels), dpi=300)
    fig.patch.set_facecolor("white")
    gs_fig = gridspec.GridSpec(n_panels, 1, hspace=0.55)

    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["font.family"]     = "sans-serif"

    fig.suptitle(
        f"GO / KEGG / Disease Ontology Enrichment — {set_label}\n"
        f"(gseapy v{gp.__version__}, live Enrichr API, BH-FDR corrected, top {top_n} terms per library)",
        fontsize=14, fontweight="bold", y=0.995, color="#0f172a"
    )

    for panel_i, (library, lib_tag, file_tag) in enumerate(LIBRARIES):
        ax = fig.add_subplot(gs_fig[panel_i])
        df = results_dict.get(lib_tag, pd.DataFrame())

        if df.empty or "Adj_P_value" not in df.columns:
            ax.text(0.5, 0.5, f"No results — {library}", ha="center",
                    va="center", transform=ax.transAxes, fontsize=11, color="#94a3b8")
            ax.set_title(f"{lib_tag} — {library}", fontsize=11, fontweight="bold", color="#0f172a")
            ax.axis("off")
            continue

        plot_df = df.copy()
        plot_df["Adj_P_value"] = pd.to_numeric(plot_df["Adj_P_value"], errors="coerce")
        plot_df = plot_df.dropna(subset=["Adj_P_value"])
        plot_df = plot_df.sort_values("Adj_P_value").head(top_n).copy()
        plot_df["-log10_FDR"] = -np.log10(plot_df["Adj_P_value"].clip(lower=1e-300))

        # Parse overlap "X/Y" → integer X
        def parse_overlap(s):
            try:
                return int(str(s).split("/")[0])
            except Exception:
                return 1

        plot_df["Overlap_N"] = plot_df["Overlap"].apply(parse_overlap)

        # Truncate long term names
        plot_df["Term_Short"] = plot_df["Term"].apply(
            lambda t: (t[:65] + "…") if len(str(t)) > 67 else t
        )

        # Sort ascending so top term is at top of axis
        plot_df = plot_df.sort_values("-log10_FDR", ascending=True)

        color = LIB_COLORS.get(lib_tag, "#64748b")
        sizes = (plot_df["Overlap_N"] / plot_df["Overlap_N"].max()) * 250 + 30

        sc = ax.scatter(
            plot_df["-log10_FDR"],
            range(len(plot_df)),
            s=sizes,
            c=color,
            alpha=0.78,
            edgecolors="white",
            linewidths=0.6,
            zorder=3,
        )

        ax.set_yticks(range(len(plot_df)))
        ax.set_yticklabels(plot_df["Term_Short"], fontsize=8.5)
        ax.set_xlabel("-log₁₀(BH-FDR Adjusted P-value)", fontsize=10, fontweight="bold")
        ax.set_title(
            f"{lib_tag} — {library}  │  {len(df)} total terms returned, top {len(plot_df)} shown",
            fontsize=10.5, fontweight="bold", color="#0f172a", pad=8
        )
        ax.axvline(-np.log10(FDR_THRESH), color="#ef4444", linewidth=1.2,
                   linestyle="--", label=f"FDR = {FDR_THRESH}")
        ax.set_xlim(left=0)
        ax.legend(fontsize=8.5, loc="lower right")
        ax.grid(axis="x", alpha=0.3, linewidth=0.7)
        ax.set_facecolor("#f8fafc")

        # Annotate FDR values on dots
        for _, row in plot_df.iterrows():
            ax.annotate(
                f"  {row['Adj_P_value']:.2e}",
                (row["-log10_FDR"], list(plot_df.index).index(row.name)),
                fontsize=7.2, va="center", color="#334155"
            )

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n  ✓ Dot-plot saved → {out_png}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN LOOP: run enrichment for each gene set, save separate files
# ─────────────────────────────────────────────────────────────────────────────
for set_key, set_info in GENE_SETS.items():
    set_tag   = set_key          # "9genes" or "5genes"
    set_label = set_info["label"]
    gene_list = set_info["genes"]

    print(f"\n{'═'*72}")
    print(f"GENE SET : {set_label}")
    print(f"Genes    : {gene_list}")
    print(f"Count    : {len(gene_list)}")
    print(f"{'═'*72}")

    results_dict = {}  # lib_tag → DataFrame

    for library, lib_tag, file_tag in LIBRARIES:
        df = run_enrichr(gene_list, library, set_tag, lib_tag)
        results_dict[lib_tag] = df

        out_csv = f"results/{file_tag}_{set_tag}.csv"
        df.to_csv(out_csv, index=False)
        print(f"     ✓ Saved {len(df)} rows → {out_csv}")

        # Brief preview of top 5 significant terms
        if not df.empty and "Adj_P_value" in df.columns:
            sig = df[pd.to_numeric(df["Adj_P_value"], errors="coerce") < FDR_THRESH]
            print(f"     Significant terms (FDR<{FDR_THRESH}): {len(sig)} / {len(df)}")
            preview = sig.sort_values("Adj_P_value").head(5)
            if not preview.empty:
                for _, row in preview.iterrows():
                    print(f"       {row['Term'][:60]:60s}  FDR={float(row['Adj_P_value']):.3e}")

        # Pause briefly between API calls to be respectful
        time.sleep(1.5)

    # Generate dot-plot figure for this gene set
    out_png = f"plots/enrichment_{set_tag}.png"
    print(f"\n  Building dot-plot figure → {out_png}")
    make_dotplot(results_dict, set_tag, set_label, out_png, top_n=TOP_N)

print(f"\n{'═'*72}")
print("ALL ENRICHMENT ANALYSES COMPLETE")
print(f"{'═'*72}")
print("\nOutputs:")
for set_key in GENE_SETS:
    for library, lib_tag, file_tag in LIBRARIES:
        print(f"  results/{file_tag}_{set_key}.csv")
    print(f"  plots/enrichment_{set_key}.png\n")
