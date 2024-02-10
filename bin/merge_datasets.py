#!/opt/conda/bin/python

import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
import numpy as np

columns_required = {
    "sex": False,
    "batch": True,
    "cell_type": False,
    "condition": False,
    "patient": True,
    "tissue": True,
    "dataset": True,
    "transfer": True
}

parser = argparse.ArgumentParser(description="Merge datasets")
parser.add_argument("--input", help="Input file", type=str, nargs="+")
parser.add_argument("--output_batches", help="Output file, batches", type=str)
parser.add_argument("--output_integration", help="Output file containing only cells which do not require transfer learning", type=str)
parser.add_argument("--output_intersection", help="Output file containing all cells but gene intersection", type=str)
parser.add_argument("--suffix_transfer", help="Output file suffix, cells for transfer learning", type=str)
parser.add_argument("--output_counts", help="Output file, outer join of cells and genes", type=str)

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
adata_outer = ad.concat(datasets, join='outer')

# Filter genes with no counts in core atlas
gene_mask, _ = sc.pp.filter_genes(adata[~adata.obs["transfer"]], min_cells=1, inplace=False)
adata = adata[:, gene_mask]

# Filter cells with no counts
cell_mask, _ = sc.pp.filter_cells(adata, min_genes=1, inplace=False)
adata = adata[cell_mask, :]
adata_outer = adata_outer[cell_mask, :]

# Filter genes with too few occurrences in outer join
sc.pp.filter_genes(adata_outer, min_cells=0.005 * adata_outer.shape[0])

# Make sure that there are no underscores in the cell names
adata.obs_names = adata.obs_names.str.replace("_", "-")
adata.obs_names_make_unique()
adata_outer.obs_names = adata.obs_names

# Convert to CSR matrix
adata.X = csr_matrix(adata.X)
adata_outer.X = csr_matrix(adata_outer.X)

adata.obs["batch"] = adata.obs["dataset"].astype(str) + "_" + adata.obs["batch"].astype(str)
adata.obs["patient"] = adata.obs["dataset"].astype(str) + "_" + adata.obs["patient"].astype(str)

def to_Florent_case(s: str):
    corrected = s.lower().strip()

    if corrected in ["na", "nan", "null", "unknown"]:
        return "Unknown"

    corrected = s \
        .replace(" ", "_") \
        .replace("-", "_")

    corrected = "".join([c if c.isalnum() or c == "_" else "" for c in corrected])

    # Make sure there is never more than one underscore
    corrected = corrected.replace("__", "_")

    if corrected.endswith("s"):
        corrected = corrected[:-1]

    corrected = corrected.strip(" _")

    if not corrected:
        return "Unknown"

    return corrected[0].upper() + corrected[1:]

for column in columns_required.keys():
    if column == "transfer":
        continue
    # Convert first to string and then to category
    adata.obs[column] = adata.obs[column].astype(str).fillna("Unknown").apply(to_Florent_case).astype("category")

adata_outer.obs = adata.obs

with open(args.output_batches, "w") as f:
    f.write("\n".join(adata[~adata.obs["transfer"]].obs["batch"].unique()))

adata.layers["counts"] = adata.X
adata_outer.layers["counts"] = adata_outer.X

adata_transfer = adata[adata.obs["transfer"]]
for dataset in adata_transfer.obs["dataset"].unique():
    adata_transfer_dataset = adata_transfer[adata_transfer.obs["dataset"] == dataset]
    adata_transfer_dataset.write_h5ad(dataset + args.suffix_transfer)

adata_notransfer = adata[~adata.obs["transfer"]]
adata_notransfer.write_h5ad(args.output_integration)
adata.write_h5ad(args.output_intersection)
adata_outer.write_h5ad(args.output_counts)