#!/opt/conda/bin/python

import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from collections import Counter

columns_required = {
    "sex": False,
    "batch": True,
    "cell_type": False,
    "condition": False,
    "patient": True,
    "tissue": True,
    "dataset": True,
}

parser = argparse.ArgumentParser(description="Merge datasets")
parser.add_argument("--input", help="Input file", type=str, nargs="+")
parser.add_argument("--output", help="Output file", type=str)

args = parser.parse_args()

datasets = [ad.read_h5ad(f) for f in args.input]

for dataset in datasets:
    # Make sure all required columns are present
    for column, required in columns_required.items():
        if column not in dataset.obs.columns:
            if required:
                raise ValueError(
                    f"Column {column} is required but not found in {dataset}"
                )
            else:
                dataset.obs[column] = "Unknown"

    # Subset columns
    dataset.obs = dataset.obs[columns_required.keys()]

adata = ad.concat(datasets)

# Perform minimal filtering to prevent NaNs
sc.pp.filter_cells(adata, min_genes=1)

# Make sure that there are no underscores in the cell names
adata.obs_names = adata.obs_names.str.replace("_", "-")

if "mito" not in adata.var.columns:
    adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# Convert to CSR matrix
adata.X = csr_matrix(adata.X)

adata.obs["batch"] = adata.obs["dataset"].astype(str) + "_" + adata.obs["batch"].astype(str)
adata.obs["patient"] = adata.obs["dataset"].astype(str) + "_" + adata.obs["patient"].astype(str)

def to_Florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "Unknown"

    corrected = s \
        .replace(" ", "_") \
        .replace("-", "_")

    if corrected.endswith("s"):
        corrected = corrected[:-1]

    corrected = corrected.strip(" _")

    if not corrected:
        return "Unknown"

    return corrected[0].upper() + corrected[1:]

for column in columns_required.keys():
    # Convert first to string and then to category
    adata.obs[column] = adata.obs[column].astype(str).fillna("Unknown").apply(to_Florent_case).astype("category")

adata.layers["counts"] = adata.X.copy()

adata.write_h5ad(args.output)