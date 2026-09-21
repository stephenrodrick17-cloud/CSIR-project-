import os
import subprocess
import numpy as np

os.makedirs('docking/receptors', exist_ok=True)

def pqr_to_pdbqt(pqr_file, pdbqt_file):
    aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
    aromatic_atoms = {
        'PHE': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TYR': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
        'TRP': {'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
        'HIS': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'}
    }

    out_lines = []
    with open(pqr_file, 'r') as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            parts = line.split()
            record = parts[0]
            try:
                serial = int(parts[1])
                name = parts[2]
                resName = parts[3]
                if len(parts) >= 11:
                    chain = parts[4]
                    resSeq = int(parts[5])
                    x = float(parts[6])
                    y = float(parts[7])
                    z = float(parts[8])
                    charge = float(parts[9])
                elif len(parts) == 10:
                    chain = 'A'
                    resSeq = int(parts[4])
                    x = float(parts[5])
                    y = float(parts[6])
                    z = float(parts[7])
                    charge = float(parts[8])
                else:
                    continue
            except Exception:
                continue

            element = name[0] if name[0].isalpha() else name[1]
            element = element.upper()
            
            if element == 'H':
                if charge < 0.12 and not ('H' in name and resName in ('LYS', 'ARG', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'CYS', 'TRP')):
                    continue
                ad_type = 'HD'
            elif element == 'C':
                if resName in aromatic_residues and name in aromatic_atoms.get(resName, set()):
                    ad_type = 'A'
                else:
                    ad_type = 'C'
            elif element == 'N':
                ad_type = 'NA' if (resName == 'HIS' and name in ('ND1', 'NE2')) else 'N'
            elif element == 'O':
                ad_type = 'OA'
            elif element == 'S':
                ad_type = 'SA'
            else:
                ad_type = element
            
            line_pdbqt = f"{record:<6s}{serial:>5d} {name:<4s} {resName:>3s} {chain:>1s}{resSeq:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}    {charge:>6.3f} {ad_type:<2s}\n"
            out_lines.append(line_pdbqt)
            
    with open(pdbqt_file, 'w') as f:
        f.writelines(out_lines)

structures = [
    ('SERPINE2_4DY0_experimental', 'docking/structures/4DY0.pdb'),
    ('SERPINE2_AF_predicted', 'docking/structures/AF-P07093-F1.pdb'),
    ('SERPINF2_AF_predicted', 'docking/structures/AF-P08697-F1.pdb'),
    ('COL1A1_1Q7D_experimental', 'docking/structures/1Q7D.pdb'),
    ('COL15A1_3N3F_experimental', 'docking/structures/3N3F.pdb')
]

receptor_configs = {}

for name, pdb_path in structures:
    pqr_path = os.path.join('docking/receptors', f"{name}.pqr")
    pdbqt_path = os.path.join('docking/receptors', f"{name}.pdbqt")
    
    print(f"\n--- Preparing receptor: {name} ---")
    cmd = f"pdb2pqr --ff AMBER --drop-water --nodebump --noopt {pdb_path} {pqr_path}"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error in pdb2pqr for {name}: {res.stderr}")
        continue
    
    pqr_to_pdbqt(pqr_path, pdbqt_path)
    
    # Compute center and span
    coords = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                p = line.split()
                coords.append([float(p[6]), float(p[7]), float(p[8])])
    arr = np.array(coords)
    center = arr.mean(axis=0).round(2)
    span = (arr.max(axis=0) - arr.min(axis=0)).round(2)
    
    # AutoDock Vina grid box: span + 4 Angstrom padding, capped appropriately
    box_size = np.clip(span + 4.0, 20.0, 45.0).round(2)
    
    receptor_configs[name] = {
        'pdbqt_path': pdbqt_path,
        'center_x': center[0],
        'center_y': center[1],
        'center_z': center[2],
        'size_x': box_size[0],
        'size_y': box_size[1],
        'size_z': box_size[2],
        'n_atoms': len(coords)
    }
    print(f"  PDBQT Atoms: {len(coords):,}")
    print(f"  Center: {center.tolist()}")
    print(f"  Grid Size: {box_size.tolist()}")

import json
with open('docking/receptors/grid_boxes.json', 'w') as f:
    json.dump(receptor_configs, f, indent=2)

print("\nAll receptor grid boxes calculated and saved to docking/receptors/grid_boxes.json!")
