#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate STRING Protein-Protein Interaction (PPI) Network Visualization
-----------------------------------------------------------------------
Creates a publication-quality PPI network for the 15 validated pan-fibrotic genes,
distinguishing ECM vs Non-ECM proteins with custom styling, and exports
Cytoscape-compatible network edge & node files (.csv).
"""

import os
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "validation_2", "network", "cytoscape_ppi")
os.makedirs(OUT_DIR, exist_ok=True)

# 15 Validated Pan-Fibrotic Genes
validated_genes_info = {
    "COL1A1": {"is_ecm": True, "category": "Collagens"},
    "COL1A2": {"is_ecm": True, "category": "Collagens"},
    "COL3A1": {"is_ecm": True, "category": "Collagens"},
    "COL15A1": {"is_ecm": True, "category": "Collagens"},
    "AEBP1": {"is_ecm": True, "category": "ECM Glycoproteins"},
    "VWF": {"is_ecm": True, "category": "ECM Glycoproteins"},
    "SERPINE2": {"is_ecm": True, "category": "ECM Regulators"},
    "CCL19": {"is_ecm": True, "category": "Secreted Factors"},
    "INMT": {"is_ecm": False, "category": "Non-ECM"},
    "LYZ": {"is_ecm": False, "category": "Non-ECM"},
    "RNASE1": {"is_ecm": False, "category": "Non-ECM"},
    "CFB": {"is_ecm": False, "category": "Non-ECM"},
    "HDAC7": {"is_ecm": False, "category": "Non-ECM"},
    "CPE": {"is_ecm": False, "category": "Non-ECM"},
    "MME": {"is_ecm": False, "category": "Non-ECM"}
}

genes = list(validated_genes_info.keys())

# Fetch STRING API edge list
string_url = "https://string-db.org/api/tsv/network"
params = {
    "identifiers": "%0d".join(genes),
    "species": 9606,
    "required_score": 400
}

url_full = f"{string_url}?identifiers={'%0d'.join(genes)}&species=9606&required_score=400"
response = requests.get(url_full)

edges = []
if response.status_code == 200:
    lines = response.text.strip().split("\n")
    if len(lines) > 1:
        header = lines[0].split("\t")
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= 6:
                p1 = parts[2]
                p2 = parts[3]
                score = float(parts[5]) if len(parts) > 5 else 0.4
                edges.append((p1, p2, score))

# Convert to DataFrame
edges_df = pd.DataFrame(edges, columns=["preferredName_A", "preferredName_B", "score"])

# Add all nodes even if isolated
nodes_df = pd.DataFrame([
    {"gene": g, "is_ecm": info["is_ecm"], "category": info["category"]}
    for g, info in validated_genes_info.items()
])

# Save Cytoscape-compatible files
edges_df.to_csv(os.path.join(OUT_DIR, "cytoscape_edges.csv"), index=False)
nodes_df.to_csv(os.path.join(OUT_DIR, "cytoscape_nodes.csv"), index=False)

print(f"Saved Cytoscape edge table ({len(edges_df)} interactions) and node table ({len(nodes_df)} genes).")

# Build NetworkX Graph
G = nx.Graph()
for g in genes:
    G.add_node(g, is_ecm=validated_genes_info[g]["is_ecm"], category=validated_genes_info[g]["category"])

for _, row in edges_df.iterrows():
    if row["preferredName_A"] in genes and row["preferredName_B"] in genes:
        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

# High quality network rendering
fig, ax = plt.subplots(figsize=(12, 10), dpi=300)

pos = nx.spring_layout(G, k=1.2, seed=42)

# Colors
ecm_color = "#e74c3c"  # Vibrant Coral/Red for ECM
non_ecm_color = "#3498db"  # Blue for Non-ECM

node_colors = [ecm_color if G.nodes[n]["is_ecm"] else non_ecm_color for n in G.nodes()]
node_sizes = [1800 if G.nodes[n]["is_ecm"] else 1200 for n in G.nodes()]

# Draw edges with opacity based on STRING combined score
if len(G.edges()) > 0:
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    widths = [w * 4.0 for w in weights]
    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, alpha=0.6, edge_color="#7f8c8d")

# Draw nodes
nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=node_sizes, edgecolors="black", linewidths=1.5)

# Draw labels
nx.draw_networkx_labels(G, pos, ax=ax, font_size=11, font_family="sans-serif", font_weight="bold", font_color="white")

# Add Custom Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='ECM Signature Genes (Collagen, Glycoproteins, Regulators)', markerfacecolor=ecm_color, markersize=15),
    Line2D([0], [0], marker='o', color='w', label='Non-ECM Validated Genes', markerfacecolor=non_ecm_color, markersize=12),
    Line2D([0], [0], color='#7f8c8d', lw=3, label='STRING Interaction (Score >= 0.400)')
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=11, frameon=True, facecolor='white', edgecolor='gray')

plt.title("Pan-Fibrotic Validated Signature: STRING Protein Interaction Network\n(ECM vs Non-ECM Functional Binding)", fontsize=14, fontweight="bold", pad=20)
plt.axis("off")
plt.tight_layout()

output_png = os.path.join(OUT_DIR, "string_ppi_ecm_vs_nonecm_network.png")
plt.savefig(output_png, dpi=300)
plt.close()

print(f"Network plot saved to: {output_png}")
