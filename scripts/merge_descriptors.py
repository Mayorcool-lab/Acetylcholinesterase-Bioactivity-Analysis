import pandas as pd

# Use snakemake variables
fp_path = snakemake.input.fingerprints
raw_path = snakemake.input.raw_data
output_path = snakemake.output.merged

# Load data
df_fp = pd.read_csv(fp_path)
df_fp = df_fp.rename(columns={'Name': 'molecule_chembl_id'})

df_raw = pd.read_csv(raw_path)
physchem_cols = ['molecule_chembl_id', 'MW', 'LogP', 'NumHDonors', 'NumHAcceptors', 'pIC50']
df_raw = df_raw[physchem_cols]

# Merge and save
df_merged = pd.merge(df_fp, df_raw, on='molecule_chembl_id')
df_merged.to_csv(output_path, index=False)
