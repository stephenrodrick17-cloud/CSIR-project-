import os
import subprocess
import json
import time

os.makedirs('docking/output', exist_ok=True)
os.makedirs('docking/logs', exist_ok=True)

with open('docking/receptors/grid_boxes.json') as f:
    grid_boxes = json.load(f)

# Docking runs: (Run_ID, Receptor_Key, Target_Gene, Target_Protein, Structure_Type, Ligand_Name, Ligand_File, Rationale)
docking_runs = [
    (
        "SERPINE2_4DY0_Camostat",
        "SERPINE2_4DY0_experimental",
        "SERPINE2",
        "Protease Nexin-1 (Glia-derived nexin)",
        "Experimental X-ray (PDB: 4DY0, 2.35 Å)",
        "Camostat",
        "docking/ligands/Camostat.pdbqt",
        "Clinically approved serine protease inhibitor; attenuates serpin-mediated fibrogenesis (Uno et al. Hepatol Res 2008)."
    ),
    (
        "SERPINE2_4DY0_Nafamostat",
        "SERPINE2_4DY0_experimental",
        "SERPINE2",
        "Protease Nexin-1 (Glia-derived nexin)",
        "Experimental X-ray (PDB: 4DY0, 2.35 Å)",
        "Nafamostat",
        "docking/ligands/Nafamostat.pdbqt",
        "Clinically approved synthetic serine protease inhibitor; blocks serpin proteolytic cascade (Mimura et al. AJP Renal 2010)."
    ),
    (
        "SERPINE2_AF_Camostat",
        "SERPINE2_AF_predicted",
        "SERPINE2",
        "Protease Nexin-1 (Glia-derived nexin)",
        "AlphaFold Predicted Model (AF-P07093-F1)",
        "Camostat",
        "docking/ligands/Camostat.pdbqt",
        "Full-length predicted model validation for Camostat binding."
    ),
    (
        "SERPINF2_AF_Nafamostat",
        "SERPINF2_AF_predicted",
        "SERPINF2",
        "Alpha-2-antiplasmin",
        "AlphaFold Predicted Model (AF-P08697-F1)",
        "Nafamostat",
        "docking/ligands/Nafamostat.pdbqt",
        "Serpin inhibitor targeting alpha-2-antiplasmin to restore fibrinolytic capacity (Iwaki et al. Kidney Int 2007)."
    ),
    (
        "SERPINF2_AF_Camostat",
        "SERPINF2_AF_predicted",
        "SERPINF2",
        "Alpha-2-antiplasmin",
        "AlphaFold Predicted Model (AF-P08697-F1)",
        "Camostat",
        "docking/ligands/Camostat.pdbqt",
        "Serine protease inhibitor evaluation against SERPINF2 predicted structure."
    ),
    (
        "COL1A1_AF_Pirfenidone",
        "COL1A1_AF_predicted",
        "COL1A1",
        "Collagen alpha-1(I) chain",
        "AlphaFold Predicted Model (AF-P02452-F1)",
        "Pirfenidone",
        "docking/ligands/Pirfenidone.pdbqt",
        "FDA-approved anti-fibrotic; directly suppresses COL1A1 collagen deposition (Noble et al. Lancet 2011)."
    ),
    (
        "COL1A1_AF_Nintedanib",
        "COL1A1_AF_predicted",
        "COL1A1",
        "Collagen alpha-1(I) chain",
        "AlphaFold Predicted Model (AF-P02452-F1)",
        "Nintedanib",
        "docking/ligands/Nintedanib.pdbqt",
        "FDA-approved anti-fibrotic; multi-kinase and procollagen fibril assembly inhibitor (Richeldi et al. NEJM 2014)."
    ),
    (
        "COL15A1_3N3F_Pirfenidone",
        "COL15A1_3N3F_experimental",
        "COL15A1",
        "Collagen alpha-1(XV) NC1 domain",
        "Experimental X-ray (PDB: 3N3F, 2.00 Å)",
        "Pirfenidone",
        "docking/ligands/Pirfenidone.pdbqt",
        "Approved anti-fibrotic tested against collagen XV NC1 non-collagenous domain."
    ),
    (
        "COL15A1_3N3F_Nintedanib",
        "COL15A1_3N3F_experimental",
        "COL15A1",
        "Collagen alpha-1(XV) NC1 domain",
        "Experimental X-ray (PDB: 3N3F, 2.00 Å)",
        "Nintedanib",
        "docking/ligands/Nintedanib.pdbqt",
        "Approved anti-fibrotic tested against collagen XV NC1 domain."
    )
]

vina_exe = os.path.abspath("bin/vina.exe")
results_summary = []

print(f"================================================================================")
print(f"STARTING GENUINE AUTODOCK VINA DOCKING SUITE (Vina v1.2.7)")
print(f"Total complexes to simulate: {len(docking_runs)}")
print(f"================================================================================\n")

for idx, (run_id, rec_key, gene, prot_desc, struct_type, lig_name, lig_path, rationale) in enumerate(docking_runs, 1):
    box = grid_boxes[rec_key]
    rec_path = os.path.abspath(box['pdbqt_path'])
    lig_full = os.path.abspath(lig_path)
    out_pdbqt = os.path.abspath(f"docking/output/{run_id}_out.pdbqt")
    log_file = os.path.abspath(f"docking/logs/{run_id}.log")
    
    # Use grid box size 30 x 30 x 30 Angstroms centered on receptor center
    cx, cy, cz = box['center_x'], box['center_y'], box['center_z']
    
    cmd_args = [
        vina_exe,
        "--receptor", rec_path,
        "--ligand", lig_full,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", "30",
        "--size_y", "30",
        "--size_z", "30",
        "--exhaustiveness", "8",
        "--num_modes", "9",
        "--out", out_pdbqt
    ]
    
    cmd_str = f'"{vina_exe}" --receptor "{rec_path}" --ligand "{lig_full}" --center_x {cx} --center_y {cy} --center_z {cz} --size_x 30 --size_y 30 --size_z 30 --exhaustiveness 8 --num_modes 9 --out "{out_pdbqt}"'
    
    print(f"[{idx}/{len(docking_runs)}] RUNNING: {run_id}")
    print(f"  Target: {gene} ({struct_type})")
    print(f"  Ligand: {lig_name}")
    print(f"  Command: {cmd_str}")
    
    t0 = time.time()
    res = subprocess.run(cmd_args, capture_output=True, text=True)
    dt = time.time() - t0
    
    # Write exact raw log output
    full_log = f"COMMAND LINE:\n{cmd_str}\n\nSTDOUT:\n{res.stdout}\n\nSTDERR:\n{res.stderr}\n"
    with open(log_file, "w") as lf:
        lf.write(full_log)
    
    # Parse best binding mode affinity from stdout
    best_affinity = None
    lines = res.stdout.splitlines()
    in_table = False
    for l in lines:
        if "-----+------------+----------+----------" in l:
            in_table = True
            continue
        if in_table:
            parts = l.split()
            if len(parts) >= 4 and parts[0] == "1":
                try:
                    best_affinity = float(parts[1])
                except ValueError:
                    pass
                break
    
    status = "COMPLETED" if res.returncode == 0 and best_affinity is not None else "FAILED"
    print(f"  Status: {status} (Elapsed: {dt:.1f}s)")
    print(f"  Top Binding Affinity (Mode 1): {best_affinity} kcal/mol")
    print(f"  Raw Log Saved: {log_file}")
    print(f"  Docked Poses Saved: {out_pdbqt}\n")
    
    results_summary.append({
        "Run_ID": run_id,
        "Target_Gene": gene,
        "Target_Protein": prot_desc,
        "Structure_Source": struct_type,
        "Drug_Candidate": lig_name,
        "Literature_Rationale": rationale,
        "Vina_Affinity_kcal_mol": best_affinity,
        "Execution_Status": status,
        "Execution_Time_sec": round(dt, 2),
        "Raw_Log_Path": log_file,
        "Docked_Poses_PDBQT": out_pdbqt,
        "Vina_Command": cmd_str
    })

import pandas as pd
df_results = pd.DataFrame(results_summary)
out_csv = "docking/results_vina_computed.csv"
df_results.to_csv(out_csv, index=False)
print(f"================================================================================")
print(f"ALL VINA SIMULATIONS FINISHED! Summary CSV saved to: {out_csv}")
print(f"================================================================================")
