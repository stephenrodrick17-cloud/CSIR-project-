# -*- coding: utf-8 -*-
"""
Master Execution Pipeline: Discovery Core, Validation 1, and Validation 2 across 4 Human Organs.
"""
import os
import glob
import pandas as pd
import numpy as np
import discovery_config

# 1. Enforce Cohort Isolation Guard
discovery_config.verify_cohort_isolation()

base_dir = r"D:\CSIR"
results_dir = os.path.join(base_dir, "results")
os.makedirs(results_dir, exist_ok=True)

# 2. Compute Discovery 4-Organ Conserved Core DEGs
kidney_deg = pd.read_csv(os.path.join(results_dir, "kidney_DEGs.csv"))
liver_deg = pd.read_csv(os.path.join(results_dir, "liver_DEGs.csv"))
lung_deg = pd.read_csv(os.path.join(results_dir, "lung_DEGs.csv"))
skin_deg = pd.read_csv(os.path.join(results_dir, "skin_DEGs.csv"))

fc_cut = 0.585
p_cut = 0.05

k_sig = kidney_deg[(kidney_deg["adj_p_value"] < p_cut) & (kidney_deg["logFC"].abs() >= fc_cut)].set_index("gene")
liv_sig = liver_deg[(liver_deg["adj_p_value"] < p_cut) & (liver_deg["logFC"].abs() >= fc_cut)].set_index("gene")
lun_sig = lung_deg[(lung_deg["adj_p_value"] < p_cut) & (lung_deg["logFC"].abs() >= fc_cut)].set_index("gene")
s_sig = skin_deg[(skin_deg["adj_p_value"] < p_cut) & (skin_deg["logFC"].abs() >= fc_cut)].set_index("gene")

core_genes = sorted(list(set(k_sig.index) & set(liv_sig.index) & set(lun_sig.index) & set(s_sig.index)))
print("=" * 80)
print(f"DISCOVERY CORE SUMMARY: {len(core_genes)} Conserved Pan-Fibrotic DEGs across all 4 organs")
print("=" * 80)

# 3. Annotate with Human Matrisome DB
matrisome_df = pd.read_excel(os.path.join(base_dir, "ECM genes all.xlsx"), header=1)
mat_dict = {}
for _, row in matrisome_df.iterrows():
    sym = str(row["Gene Symbol"]).strip().upper()
    cat = str(row["Matrisome Category"]).strip()
    div = str(row["Matrisome Division"]).strip()
    if sym and sym not in ("NAN", "NA", "-"):
        mat_dict[sym] = {"category": cat, "division": div}

core_records = []
for g in core_genes:
    is_ecm = g in mat_dict
    cat = mat_dict[g]["category"] if is_ecm else "Non-ECM"
    div = mat_dict[g]["division"] if is_ecm else "Non-ECM"
    core_records.append({
        "gene": g,
        "kidney_logFC": k_sig.loc[g, "logFC"],
        "kidney_adj_p_value": k_sig.loc[g, "adj_p_value"],
        "liver_logFC": liv_sig.loc[g, "logFC"],
        "liver_adj_p_value": liv_sig.loc[g, "adj_p_value"],
        "lung_logFC": lun_sig.loc[g, "logFC"],
        "lung_adj_p_value": lun_sig.loc[g, "adj_p_value"],
        "skin_logFC": s_sig.loc[g, "logFC"],
        "skin_adj_p_value": s_sig.loc[g, "adj_p_value"],
        "is_ECM_gene": is_ecm,
        "matrisome_category": cat,
        "matrisome_division": div,
    })

df_core = pd.DataFrame(core_records)
df_core.to_csv(os.path.join(results_dir, "pan_fibrotic_core_genes_corrected.csv"), index=False)

df_ecm = df_core[df_core["is_ECM_gene"] == True].reset_index(drop=True)
df_ecm.to_csv(os.path.join(results_dir, "ecm_clean_genes.csv"), index=False)

print(f"-> Saved pan_fibrotic_core_genes_corrected.csv ({len(df_core)} genes)")
print(f"-> Saved ecm_clean_genes.csv ({len(df_ecm)} clean core ECM genes)")
print(f"Clean ECM Genes ({len(df_ecm)}): {df_ecm['gene'].tolist()}\n")

# Global bioDBnet mapping for probe conversion
global_bio_map = discovery_config.verify_cohort_isolation and {}
for f in glob.glob(os.path.join(base_dir, "**", "bioDBnet*.txt"), recursive=True):
    with open(f, "r", encoding="utf-8", errors="ignore") as fp:
        for line in fp:
            parts = line.strip().split("\t")
            if len(parts) >= 2 and parts[1].strip().upper() not in ("-", "NONE", "NAN", "NA"):
                global_bio_map[parts[0].strip()] = parts[1].strip().upper()

acc_map_file = os.path.join(base_dir, "geo_cache", "gb_acc_to_symbol_map.json")
if os.path.exists(acc_map_file):
    import json
    with open(acc_map_file, "r") as f:
        global_bio_map.update(json.load(f))


def prep_table(fpath, is_val1_liver=False):
    df = pd.read_csv(fpath, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    cols_l = {c.lower(): c for c in df.columns}
    
    lfc_col = next((cols_l[c] for c in ["logfc", "log2foldchange", "lfc"] if c in cols_l), None)
    p_col = next((cols_l[c] for c in ["adj.p.val", "padj", "adj_p_value", "fdr"] if c in cols_l), None)
    sym_col = next((cols_l[c] for c in ["symbol", "gene.symbol", "gene_symbol", "genesymbol"] if c in cols_l), None)
    id_col = cols_l.get("id", df.columns[0])
    
    if sym_col:
        df["gene"] = df[sym_col].astype(str).str.strip().str.upper()
    else:
        df["gene"] = df[id_col].astype(str).str.strip().map(global_bio_map)
        df["gene"] = df["gene"].fillna(df[id_col].astype(str).str.strip().str.upper())
        
    df["lfc_clean"] = pd.to_numeric(df[lfc_col], errors="coerce")
    if is_val1_liver:
        df["lfc_clean"] = -df["lfc_clean"]  # Invert Control vs Cirrhotic
        
    df["padj_clean"] = pd.to_numeric(df[p_col], errors="coerce")
    return df.dropna(subset=["gene", "lfc_clean", "padj_clean"]).sort_values("padj_clean").drop_duplicates("gene", keep="first").set_index("gene")


# 4. Process Validation 1 Datasets
print("=" * 80)
print("VALIDATION LAYER 1: 4-ORGAN INDEPENDENT COHORTS")
print("=" * 80)
val1_k = prep_table(r"D:\CSIR\Kidney\18\GSE200818.top.table.tsv")
val1_l = prep_table(r"D:\CSIR\Liver\694\GSE162694.top.table.tsv", is_val1_liver=True)
val1_lu = prep_table(r"D:\CSIR\Lungs\206\GSE24206.top.table (2).tsv")
val1_s = prep_table(r"D:\CSIR\Skin\095\GSE58095.top.table.tsv")

val1_records = []
for g in df_ecm["gene"]:
    k_sig = (val1_k.loc[g, "padj_clean"] < 0.05) if g in val1_k.index else False
    l_sig = (val1_l.loc[g, "padj_clean"] < 0.05) if g in val1_l.index else False
    lu_sig = (val1_lu.loc[g, "padj_clean"] < 0.05) if g in val1_lu.index else False
    s_sig = (val1_s.loc[g, "padj_clean"] < 0.05) if g in val1_s.index else False
    
    n_sig = sum([k_sig, l_sig, lu_sig, s_sig])
    val1_records.append({
        "gene": g,
        "kidney_val1_sig": k_sig, "liver_val1_sig": l_sig,
        "lung_val1_sig": lu_sig, "skin_val1_sig": s_sig,
        "val1_n_sig_organs": n_sig
    })

df_val1_res = pd.DataFrame(val1_records)
df_val1_res.to_csv(os.path.join(results_dir, "validation1_layer2_all_results.csv"), index=False)
n1_ge1 = (df_val1_res["val1_n_sig_organs"] >= 1).sum()
n1_ge2 = (df_val1_res["val1_n_sig_organs"] >= 2).sum()
print(f"Val 1 Replication in >= 1 Organ: {n1_ge1}/{len(df_val1_res)} ({n1_ge1/len(df_val1_res)*100:.1f}%)")
print(f"Val 1 Replication in >= 2 Organs: {n1_ge2}/{len(df_val1_res)} ({n1_ge2/len(df_val1_res)*100:.1f}%)\n")

# 5. Process Validation 2 Datasets
print("=" * 80)
print("VALIDATION LAYER 2: 4-ORGAN HELD-OUT MULTI-PLATFORM COHORTS")
print("=" * 80)
val2_k = prep_table(r"D:\CSIR\Kidney\Validation 2\GSE30529.top.table.tsv")
val2_l = prep_table(r"D:\CSIR\Liver\validate 2\GSE14323.top.table.tsv")
val2_lu = prep_table(r"D:\CSIR\Lungs\validate 2\GSE83717.top.table.tsv")
val2_s = prep_table(r"D:\CSIR\Skin\validation 2\GSE125362.top.table.tsv")

val2_records = []
for g in df_ecm["gene"]:
    k_sig = (val2_k.loc[g, "padj_clean"] < 0.05) if g in val2_k.index else False
    l_sig = (val2_l.loc[g, "padj_clean"] < 0.05) if g in val2_l.index else False
    lu_sig = (val2_lu.loc[g, "padj_clean"] < 0.05) if g in val2_lu.index else False
    s_sig = (val2_s.loc[g, "padj_clean"] < 0.05) if g in val2_s.index else False
    
    n_sig = sum([k_sig, l_sig, lu_sig, s_sig])
    val2_records.append({
        "gene": g,
        "kidney_val2_sig": k_sig, "liver_val2_sig": l_sig,
        "lung_val2_sig": lu_sig, "skin_val2_sig": s_sig,
        "val2_n_sig_organs": n_sig
    })

df_val2_res = pd.DataFrame(val2_records)
df_val2_res.to_csv(os.path.join(results_dir, "validation2_ecm_core_results.csv"), index=False)
n2_ge1 = (df_val2_res["val2_n_sig_organs"] >= 1).sum()
n2_ge2 = (df_val2_res["val2_n_sig_organs"] >= 2).sum()
print(f"Val 2 Replication in >= 1 Organ: {n2_ge1}/{len(df_val2_res)} ({n2_ge1/len(df_val2_res)*100:.1f}%)")
print(f"Val 2 Replication in >= 2 Organs: {n2_ge2}/{len(df_val2_res)} ({n2_ge2/len(df_val2_res)*100:.1f}%)")

# Detailed per-gene summary table
df_summary_all = pd.merge(df_val1_res[["gene", "val1_n_sig_organs"]], df_val2_res, on="gene")
print("\n" + "=" * 80)
print("FINAL REPLICATION SUMMARY FOR 24 CORE ECM GENES")
print("=" * 80)
print(df_summary_all.sort_values(by=["val2_n_sig_organs", "val1_n_sig_organs"], ascending=False).to_string(index=False))