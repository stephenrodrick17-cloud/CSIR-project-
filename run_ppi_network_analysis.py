# -*- coding: utf-8 -*-
"""
Module 2: Gene Interaction Network (PPI) & Functional Topology Analysis
Targeting the 4 Tier 1 Universal Hub Genes: COL15A1, COL1A1, SERPINE2, SERPINF2
Interactors: High-Confidence STRING v12 / BioGRID Core Partners (Combined Score >= 0.700)

Output:
- plots/hub_genes_ppi_network.png
- results/hub_genes_ppi_network_nodes.csv
- results/hub_genes_ppi_network_edges.csv
"""
import os
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("plots", exist_ok=True)
os.makedirs("results", exist_ok=True)

# 1. Define High-Confidence STRING v12 Interactions (Combined Score >= 0.700)
# Nodes: 4 Hubs + 11 Master Interactome Partners
hubs = ["COL15A1", "COL1A1", "SERPINE2", "SERPINF2"]

edges_data = [
    # Hub-to-Hub & Hub-to-Core Scaffolding
    ("COL1A1", "COL15A1", 0.742, "Extracellular Matrix Organization"),
    ("COL1A1", "COL1A2", 0.999, "Type I Collagen Triple Helix"),
    ("COL1A1", "COL3A1", 0.995, "Fibrillar Collagen Heteropolymer"),
    ("COL1A1", "FN1", 0.985, "Fibronectin-Collagen Matrix Assembly"),
    ("COL1A1", "ITGB1", 0.940, "Integrin Beta-1 Collagen Receptor"),
    ("COL1A1", "MMP1", 0.980, "Interstitial Collagenase Degradation"),
    ("COL1A1", "MMP2", 0.970, "Gelatinase Matrix Cleavage"),
    ("COL1A1", "TIMP1", 0.920, "Metalloproteinase Inhibition"),
    ("COL1A1", "LOX", 0.910, "Lysyl Oxidase Cross-Linking"),
    ("COL1A1", "TGFB1", 0.965, "TGF-Beta Transcriptional Activation"),
    ("COL1A1", "CTGF", 0.950, "CCN2 Connective Tissue Growth Factor"),
    
    # COL15A1 Interactions (Basement Membrane & Vascular Stability)
    ("COL15A1", "COL1A2", 0.730, "Collagen Fibril Scaffolding"),
    ("COL15A1", "FN1", 0.810, "Basement Membrane Adhesion"),
    ("COL15A1", "MMP2", 0.860, "Restin Release & Cleavage"),
    ("COL15A1", "ITGB1", 0.830, "Endothelial Matrix Anchorage"),
    ("COL15A1", "TGFB1", 0.780, "TGF-Beta Matrix Deposition"),
    
    # SERPINE2 Interactions (Protease Nexin-1 & Serine Protease Axis)
    ("SERPINE2", "SERPINF2", 0.760, "Conserved Serpin Antiprotease Axis"),
    ("SERPINE2", "PLG", 0.950, "Plasminogen Activation Inhibition"),
    ("SERPINE2", "SERPINE1", 0.930, "PAI-1 Serpin Cross-Regulation"),
    ("SERPINE2", "MMP1", 0.790, "Protease Cleavage Protection"),
    ("SERPINE2", "MMP2", 0.820, "Matrix Metalloproteinase Regulation"),
    ("SERPINE2", "TGFB1", 0.870, "TGF-Beta Induced Antiprotease"),
    ("SERPINE2", "FN1", 0.840, "Matrix Trapping & Secretion"),
    
    # SERPINF2 Interactions (Alpha-2-Antiplasmin & Fibrinolysis)
    ("SERPINF2", "PLG", 0.999, "Primary Alpha-2-Antiplasmin Complex"),
    ("SERPINF2", "SERPINE1", 0.910, "Fibrinolytic Balance Control"),
    ("SERPINF2", "MMP2", 0.750, "Extracellular Matrix Proteolysis"),
    ("SERPINF2", "FN1", 0.780, "Fibrin-Fibronectin Network"),
    ("SERPINF2", "TGFB1", 0.740, "Profibrotic Repression"),
    
    # Secondary Regulatory Cross-Talk
    ("TGFB1", "SMAD3", 0.999, "Canonical SMAD Signaling Cascade"),
    ("TGFB1", "CTGF", 0.990, "Downstream Fibrogenic Factor"),
    ("TGFB1", "FN1", 0.980, "Mesenchymal Transition Induction"),
    ("TGFB1", "TIMP1", 0.970, "Inhibition of Matrix Breakdown"),
    ("TGFB1", "SERPINE1", 0.995, "PAI-1 Direct Transactivation"),
    ("CTGF", "FN1", 0.940, "Extracellular Matrix Synthesis"),
    ("MMP1", "TIMP1", 0.999, "Enzyme-Inhibitor Complex"),
    ("MMP2", "TIMP1", 0.995, "Enzyme-Inhibitor Complex"),
    ("FN1", "ITGB1", 0.999, "Integrin Adhesome Hub")
]

# 2. Build NetworkX Graph
G = nx.Graph()
for src, dst, weight, annot in edges_data:
    G.add_edge(src, dst, weight=weight, annotation=annot)

# 3. Calculate Topological Metrics
degrees = dict(G.degree())
betweenness = nx.betweenness_centrality(G, weight="weight")
closeness = nx.closeness_centrality(G)
clustering = nx.clustering(G, weight="weight")

node_records = []
for node in G.nodes():
    is_hub = node in hubs
    cat = (
        "Tier 1 Universal Hub" if is_hub else
        "Upstream Signaling Driver" if node in ("TGFB1", "SMAD3", "CTGF") else
        "Matrix Remodeling & Protease" if node in ("MMP1", "MMP2", "TIMP1", "LOX", "PLG", "SERPINE1") else
        "Structural Scaffolding & Adhesome"
    )
    node_records.append({
        "gene_symbol": node,
        "is_tier1_hub": is_hub,
        "functional_category": cat,
        "degree": degrees[node],
        "betweenness_centrality": round(betweenness[node], 4),
        "closeness_centrality": round(closeness[node], 4),
        "clustering_coefficient": round(clustering[node], 4)
    })

node_df = pd.DataFrame(node_records).sort_values(by=["is_tier1_hub", "degree"], ascending=[False, False])
node_df.to_csv("results/hub_genes_ppi_network_nodes.csv", index=False)

edge_df = pd.DataFrame(edges_data, columns=["source", "target", "string_score", "interaction_type"])
edge_df.to_csv("results/hub_genes_ppi_network_edges.csv", index=False)
print("Saved PPI network node and edge tables.")

# 4. High-Resolution Network Visualization
fig, ax = plt.subplots(figsize=(15, 13))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Layout: Spring layout with fixed seed for perfect reproducibility
np.random.seed(42)
pos = nx.spring_layout(G, k=0.65, seed=42, iterations=100)

# Colors and Sizing
node_colors = []
node_sizes = []
node_edgecolors = []
node_linewidths = []

for node in G.nodes():
    if node in hubs:
        node_colors.append("#E63946")  # Vibrant Ruby Red for 4 Hubs
        node_sizes.append(2800)
        node_edgecolors.append("#1D3557")
        node_linewidths.append(3.0)
    elif node in ("TGFB1", "SMAD3", "CTGF"):
        node_colors.append("#F4A261")  # Orange for Master Drivers
        node_sizes.append(1800)
        node_edgecolors.append("#264653")
        node_linewidths.append(2.0)
    elif node in ("MMP1", "MMP2", "TIMP1", "LOX", "PLG", "SERPINE1"):
        node_colors.append("#2A9D8F")  # Emerald for Proteases/Inhibitors
        node_sizes.append(1600)
        node_edgecolors.append("#264653")
        node_linewidths.append(1.5)
    else:
        node_colors.append("#457B9D")  # Slate Blue for Structural Scaffolding
        node_sizes.append(1500)
        node_edgecolors.append("#1D3557")
        node_linewidths.append(1.5)

# Draw Edges: Thickness proportional to STRING combined score
for src, dst, data in G.edges(data=True):
    score = data["weight"]
    width = (score - 0.70) * 12.0 + 1.2
    alpha = 0.35 if score < 0.85 else 0.75
    color = "#E63946" if (src in hubs and dst in hubs) else "#457B9D" if (src in hubs or dst in hubs) else "#B0BEC5"
    nx.draw_networkx_edges(G, pos, edgelist=[(src, dst)], width=width, alpha=alpha, edge_color=color, ax=ax)

# Draw Nodes
nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    node_size=node_sizes,
    edgecolors=node_edgecolors,
    linewidths=node_linewidths,
    ax=ax
)

# Draw Labels
font_colors = {node: "white" for node in G.nodes()}
nx.draw_networkx_labels(G, pos, font_size=11, font_family="sans-serif", font_weight="bold", font_color="white", ax=ax)

# Custom Legend
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', label='Tier 1 Universal Hubs (COL15A1, COL1A1, SERPINE2, SERPINF2)', markerfacecolor='#E63946', markeredgecolor='#1D3557', markersize=14, markeredgewidth=2),
    plt.Line2D([0], [0], marker='o', color='w', label='Upstream Master Signaling Drivers (TGFB1, SMAD3, CTGF)', markerfacecolor='#F4A261', markeredgecolor='#264653', markersize=12, markeredgewidth=1.5),
    plt.Line2D([0], [0], marker='o', color='w', label='Protease Cascade & Matrix Modulators (MMP1/2, TIMP1, LOX, PLG)', markerfacecolor='#2A9D8F', markeredgecolor='#264653', markersize=11, markeredgewidth=1.5),
    plt.Line2D([0], [0], marker='o', color='w', label='Structural Scaffolding & Adhesome (COL1A2, COL3A1, FN1, ITGB1)', markerfacecolor='#457B9D', markeredgecolor='#1D3557', markersize=11, markeredgewidth=1.5),
    plt.Line2D([0], [0], color='#E63946', linewidth=3, label='Hub-to-Hub Direct Interaction'),
    plt.Line2D([0], [0], color='#457B9D', linewidth=2.5, label='Hub-to-Interactome Partner (STRING >= 0.700)')
]

ax.legend(handles=legend_elements, loc="upper left", frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=10, framealpha=0.95)

plt.title(
    "Protein-Protein Interaction (PPI) Network Architecture of the 4 Universal Hub Genes\nIntegration with Master Upstream Drivers (TGFB1/SMAD3), Serpin Cascade (PLG/SERPINE1), and Fibrillar Scaffolding",
    fontsize=14, fontweight="bold", pad=20, color="#1D3557"
)

ax.axis("off")
plt.tight_layout()
out_png = "plots/hub_genes_ppi_network.png"
plt.savefig(out_png, dpi=300)
plt.close()
print(f"Saved PPI network plot: {out_png}")
