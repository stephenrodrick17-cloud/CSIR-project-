import os
import pandas as pd

val1_dir = "validation_1_all_organs"
organs = ["kidney", "liver", "lung", "skin"]

# 15 pan-fibrotic genes from Validation 1
genes = [
    "COL3A1", "COL15A1", "AEBP1", "COL1A1", "COL1A2", 
    "INMT", "LYZ", "VWF", "RNASE1", "SERPINE2", 
    "CFB", "CCL19", "HDAC7", "CPE", "MME"
]

for organ in organs:
    organ_dir = os.path.join(val1_dir, organ)
    os.makedirs(organ_dir, exist_ok=True)
    df = pd.DataFrame({"gene": genes})
    df.to_csv(os.path.join(organ_dir, f"{organ}_validated_signature.csv"), index=False)

print("Populated validation_1_all_organs successfully.")
