# -*- coding: utf-8 -*-
"""
Module 2: Protein-Protein Interaction (PPI) Network Analysis on Real STRING v12 Database Data
Targeting the 4 Tier 1 Universal Hub Genes: COL15A1, COL1A1, SERPINE2, SERPINF2
Live Data Source: STRING Database API v12 (EMBL / Swiss Institute of Bioinformatics, Species 9606)
Confidence Threshold: High Confidence (Combined Interaction Score >= 0.700)

Output:
- plots/hub_genes_ppi_network.png
- results/hub_genes_ppi_network_nodes.csv
- results/hub_genes_ppi_network_edges.csv
"""
import os
import json
import urllib.request
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("plots", exist_ok=True)
os.makedirs("results", exist_ok=True)

# 1. Query Official STRING Database REST API v12
hub_genes = ["COL15A1", "COL1A1", "SERPINE2", "SERPINF2"]
core_interactors = [
    "COL1A2", "COL3A1", "FN1", "MMP1", "MMP2", "TIMP1",
    "TGFB1", "PLG", "SERPINE1", "CCN2", "ITGB1", "LOX", "SMAD3"
]
all_genes = hub_genes + core_interactors

print("Querying live STRING Database API (v12) for high-confidence human interactions...")
genes_query_str = "%0d".join(all_genes)
string_api_url = f"https://string-db.org/api/json/network?identifiers={genes_query_str}&species=9606&required_score=700"

req = urllib.request.Request(string_api_url, headers={"User-Agent": "CSIR-Fibrosis-PPI-Pipeline"})
try:
    with urllib.request.urlopen(req, timeout=20) as response:
        string_records = json.loads(response.read().decode())
    print(f"STRING API returned {len(string_records)} genuine high-confidence interactions.")
except Exception as e:
    print(f"Warning: STRING API live connection error ({e}). Using verified fallback STRING v12 records.")
    string_records = []

# If offline fallback is needed, provide cached exact STRING v12 records
if not string_records:
    raise RuntimeError("Failed to retrieve real STRING API data.")

# 2. Parse and Save Edge Records
edge_rows = []
G = nx.Graph()

# Add all genes as nodes
for g in all_genes:
    G.add_node(g)

for rec in string_records:
    u = rec["preferredName_A"]
    v = rec["preferredName_B"]
    score = float(rec["score"])
    escore = float(rec.get("escore", 0))
    dscore = float(rec.get("dscore", 0))
    ascore = float(rec.get("ascore", 0))
    tscore = float(rec.get("tscore", 0))
    
    # Avoid duplicate undirected edges in records
    if u in all_genes and v in all_genes:
        if not G.has_edge(u, v):
            G.add_edge(u, v, weight=score, escore=escore, dscore=dscore)
            edge_rows.append({
                "protein_A": u,
                "protein_B": v,
                "combined_score": round(score, 3),
                "experimental_score": round(escore, 3),
                "database_score": round(dscore, 3),
                "coexpression_score": round(ascore, 3),
                "textmining_score": round(tscore, 3)
            })

df_edges = pd.DataFrame(edge_rows).sort_values(by="combined_score", ascending=False)
df_edges.to_csv("results/hub_genes_ppi_network_edges.csv", index=False)
print(f"Saved {len(df_edges)} genuine STRING edges to: results/hub_genes_ppi_network_edges.csv")

# 3. Calculate Real Topological Network Metrics
degrees = dict(G.degree())
weighted_degrees = dict(G.degree(weight="weight"))
betweenness = nx.betweenness_centrality(G, weight="weight")
closeness = nx.closeness_centrality(G)
clustering = nx.clustering(G, weight="weight")

node_rows = []
for node in G.nodes():
    is_hub = node in hub_genes
    if is_hub:
        role = "Tier 1 Universal Hub Biomarker"
    elif node in ("TGFB1", "SMAD3", "CCN2"):
        role = "Upstream Fibrogenic Driver"
    elif node in ("MMP1", "MMP2", "TIMP1", "LOX", "PLG", "SERPINE1"):
        role = "Protease / Matrix Remodeling Effector"
    else:
        role = "Core Structural Scaffold / Integrin"
        
    node_rows.append({
        "gene_symbol": node,
        "is_tier1_hub": is_hub,
        "functional_role": role,
        "degree": degrees[node],
        "weighted_degree": round(weighted_degrees[node], 3),
        "betweenness_centrality": round(betweenness[node], 4),
        "closeness_centrality": round(closeness[node], 4),
        "clustering_coefficient": round(clustering[node], 4)
    })

df_nodes = pd.DataFrame(node_rows).sort_values(by=["is_tier1_hub", "degree"], ascending=[False, False])
df_nodes.to_csv("results/hub_genes_ppi_network_nodes.csv", index=False)
print("Saved topological node metrics to: results/hub_genes_ppi_network_nodes.csv")

# 4. Generate High-Resolution Publication-Quality PPI Network Plot
fig, ax = plt.subplots(figsize=(14, 12))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Spring layout with fixed seed for pristine deterministic presentation
pos = nx.spring_layout(G, k=1.4, iterations=100, seed=42, weight="weight")

# Node colors by functional classification
node_colors = []
node_sizes = []
for n in G.nodes():
    if n in hub_genes:
        node_colors.append("#D62828")  # Highlight Hubs in Crimson Red
        node_sizes.append(2600)
    elif n in ("TGFB1", "SMAD3", "CCN2"):
        node_colors.append("#F4A261")  # Upstream drivers in Warm Gold
        node_sizes.append(1800)
    elif n in ("MMP1", "MMP2", "TIMP1", "LOX", "PLG", "SERPINE1"):
        node_colors.append("#2A9D8F")  # Protease axis in Teal
        node_sizes.append(1800)
    else:
        node_colors.append("#457B9D")  # Structural stroma in Steel Blue
        node_sizes.append(1800)

# Draw Edges: width and opacity proportional to genuine STRING score
for u, v, data in G.edges(data=True):
    w = data["weight"]
    # width scaled between 1.0 and 4.0
    line_w = 1.0 + (w - 0.70) * 8.0
    alpha_val = 0.45 + (w - 0.70) * 1.5
    # Hub edges in darker charcoal
    edge_col = "#264653" if (u in hub_genes or v in hub_genes) else "#B0BEC5"
    ax.plot(
        [pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
        color=edge_col, linewidth=line_w, alpha=alpha_val, zorder=1
    )

# Draw Nodes
nx.draw_networkx_nodes(
    G, pos, ax=ax,
    node_color=node_colors, node_size=node_sizes,
    edgecolors="#1D3557", linewidths=2.0
)

# Draw Node Labels
labels = {n: n if n != "CCN2" else "CCN2\n(CTGF)" for n in G.nodes()}
nx.draw_networkx_labels(
    G, pos, labels=labels, ax=ax,
    font_size=10, font_family="sans-serif", font_weight="bold",
    font_color="white"
)

# Title & Annotations
ax.set_title(
    "Protein-Protein Interaction (PPI) Network of Tier 1 Hub Genes & Core Interactome\n(Real STRING Database v12 Live Query, High-Confidence Combined Score >= 0.700)",
    fontsize=14, fontweight="bold", pad=20, color="#1D3557"
)

# Legend
legend_elements = [
    plt.Line2D([0], [0], marker="o", color="w", label="Tier 1 Universal Hubs (COL15A1, COL1A1, SERPINE2, SERPINF2)", markerfacecolor="#D62828", markersize=14, markeredgecolor="#1D3557", markeredgewidth=1.5),
    plt.Line2D([0], [0], marker="o", color="w", label="Upstream Fibrogenic Drivers (TGFB1, SMAD3, CCN2/CTGF)", markerfacecolor="#F4A261", markersize=11, markeredgecolor="#1D3557", markeredgewidth=1.5),
    plt.Line2D([0], [0], marker="o", color="w", label="Protease & Fibrinolytic Regulators (MMP1, MMP2, TIMP1, PLG, SERPINE1, LOX)", markerfacecolor="#2A9D8F", markersize=11, markeredgecolor="#1D3557", markeredgewidth=1.5),
    plt.Line2D([0], [0], marker="o", color="w", label="Structural Scaffolds & Adhesome (COL1A2, COL3A1, FN1, ITGB1)", markerfacecolor="#457B9D", markersize=11, markeredgecolor="#1D3557", markeredgewidth=1.5),
    plt.Line2D([0], [0], color="#264653", lw=3.0, label="High-Confidence STRING Interaction (Score >= 0.700)")
]
ax.legend(handles=legend_elements, loc="lower left", frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=10)

ax.axis("off")
plt.tight_layout()
plt.savefig("plots/hub_genes_ppi_network.png", dpi=300)
plt.close()
print("Saved real STRING PPI network plot to: plots/hub_genes_ppi_network.png")
