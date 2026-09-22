# -*- coding: utf-8 -*-
"""
Module 2: Protein-Protein Interaction (PPI) Network Analysis on Real STRING v12 Database Data
Targeting the 5 Validated Pan-Fibrotic Hub Genes: COL15A1, COL1A1, COL3A1, SERPINE2, SERPINF2
Live Data Source: STRING Database API v12 (Species 9606, Homo sapiens)
Confidence Threshold: High Confidence (Combined Interaction Score >= 0.700)

Output:
- plots/hub_genes_ppi_network_5genes.png
- results/hub_genes_ppi_network_nodes_5genes.csv
- results/hub_genes_ppi_network_edges_5genes.csv
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
hub_genes = ["COL15A1", "COL1A1", "COL3A1", "SERPINE2", "SERPINF2"]
core_interactors = [
    "COL1A2", "FN1", "MMP1", "MMP2", "TIMP1",
    "TGFB1", "PLG", "SERPINE1", "CCN2", "ITGB1", "LOX", "SMAD3"
]
all_genes = hub_genes + core_interactors

print("Querying live STRING Database API (v12) for 5 hub genes...")
genes_query_str = "%0d".join(all_genes)
string_api_url = f"https://string-db.org/api/json/network?identifiers={genes_query_str}&species=9606&required_score=700"

req = urllib.request.Request(string_api_url, headers={"User-Agent": "CSIR-Fibrosis-PPI-Pipeline"})
try:
    with urllib.request.urlopen(req, timeout=20) as response:
        string_records = json.loads(response.read().decode())
    print(f"STRING API returned {len(string_records)} genuine high-confidence interactions.")
except Exception as e:
    print(f"Error connecting to STRING API: {e}")
    raise RuntimeError("Failed to retrieve real STRING API data.")

# 2. Parse and Save Edge Records
edge_rows = []
G = nx.Graph()

for g in all_genes:
    G.add_node(g)

for rec in string_records:
    u = rec["preferredName_A"]
    v = rec["preferredName_B"]
    score = float(rec["score"])
    
    if u in all_genes and v in all_genes and u != v:
        G.add_edge(u, v, weight=score)
        edge_rows.append({
            "source": u,
            "target": v,
            "combined_score": score,
            "interaction_type": "High-Confidence STRING Physical/Functional Association"
        })

df_edges = pd.DataFrame(edge_rows).drop_duplicates(subset=["source", "target"])
df_edges.to_csv("results/hub_genes_ppi_network_edges_5genes.csv", index=False)
print(f"Saved: results/hub_genes_ppi_network_edges_5genes.csv ({len(df_edges)} unique edges)")

# 3. Calculate Real Topological Centralities
degrees = dict(G.degree())
betweenness = nx.betweenness_centrality(G)
closeness = nx.closeness_centrality(G)

node_rows = []
for n in G.nodes():
    node_rows.append({
        "gene_symbol": n,
        "role": "Validated Hub Gene" if n in hub_genes else "Functional Partner / Upstream Master Regulator",
        "degree": degrees[n],
        "betweenness_centrality": round(betweenness[n], 4),
        "closeness_centrality": round(closeness[n], 4)
    })

df_nodes = pd.DataFrame(node_rows).sort_values(by="degree", ascending=False)
df_nodes.to_csv("results/hub_genes_ppi_network_nodes_5genes.csv", index=False)
print("Saved: results/hub_genes_ppi_network_nodes_5genes.csv")

# 4. Generate Publication-Quality Network Figure
fig, ax = plt.subplots(figsize=(13, 11), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

pos = nx.spring_layout(G, seed=42, k=0.72, iterations=120)

# Colors and Sizes
node_colors = []
node_sizes = []
for n in G.nodes():
    if n in hub_genes:
        node_colors.append("#dc2626")  # Crimson Red for Hub Genes
        node_sizes.append(1850)
    else:
        node_colors.append("#3b82f6")  # Blue for Core Interactors
        node_sizes.append(1150)

# Draw Edges
weights = [G[u][v]["weight"] * 2.2 for u, v in G.edges()]
nx.draw_networkx_edges(G, pos, ax=ax, width=weights, alpha=0.35, edge_color="#64748b")

# Draw Nodes
nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=node_sizes, edgecolors="#0f172a", linewidths=1.5)

# Labels
labels = {n: n for n in G.nodes()}
nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=7.2, font_family="sans-serif", font_weight="bold", font_color="#ffffff")

# Custom Legend
hub_patch = plt.Line2D([0], [0], marker='o', color='w', label='Validated Hub Gene (5 Consensus Biomarkers)',
                       markerfacecolor='#dc2626', markersize=14, markeredgecolor='#0f172a', markeredgewidth=1.2)
partner_patch = plt.Line2D([0], [0], marker='o', color='w', label='Interacting Driver (Proteases / Signaling / Scaffolding)',
                           markerfacecolor='#3b82f6', markersize=11, markeredgecolor='#0f172a', markeredgewidth=1.2)
ax.legend(handles=[hub_patch, partner_patch], loc="lower left", frameon=True, facecolor="white", edgecolor="#cbd5e1", fontsize=9.5)

ax.set_title("Protein-Protein Interaction (PPI) Network Architecture (5 Validated Hubs)\nLive STRING Database v12 High-Confidence Interactome (Combined Score ≥ 0.700)",
             fontsize=13.5, fontweight="bold", pad=15, color="#0f172a")

ax.axis("off")
plt.tight_layout()
out_png = "plots/hub_genes_ppi_network_5genes.png"
plt.savefig(out_png, dpi=300)
plt.close()
print(f"Saved 5-hub PPI Network figure: {out_png}")
