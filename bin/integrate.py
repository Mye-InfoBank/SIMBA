#!/usr/bin/env python3

import anndata as ad
import pandas as pd
import scib
import scanpy as sc
import argparse
from threadpoolctl import threadpool_limits
import torch

# Disable warnings
import warnings
warnings.filterwarnings("ignore")

torch.set_float32_matmul_precision('medium')

methods = {
    "bbknn": scib.ig.bbknn,
    "combat": scib.ig.combat,
    "desc": scib.ig.desc,
    "harmony": scib.ig.harmony,
    "mnn": scib.ig.mnn,
    #"saucie": scib.ig.saucie, # Deactivated because dependencies are not available
    "scanorama": scib.ig.scanorama,
    "scanvi": scib.ig.scanvi,
    "scgen": scib.ig.scgen,
    "scvi": scib.ig.scvi,
    #"trvae": scib.ig.trvae, # Deactivated because it requires keras==2.2.4 which is incompatible with other dependencies
    "trvaep": scib.ig.trvaep,
}

parser = argparse.ArgumentParser(description='Integrate data')
parser.add_argument('--input', help='Input file', type=str)
parser.add_argument('--output', help='Output file', type=str)
parser.add_argument('--hvgs', help='Pickle containing highly variable genes', type=str)
parser.add_argument('--method', help='Integration method', type=str, choices=list(methods.keys()) + ["unintegrated"])
parser.add_argument('--scvi_model', help='scvi model', type=str)
parser.add_argument('--cpus', help='Number of cpus', type=int, default=1)

args = parser.parse_args()

threadpool_limits(args.cpus)
sc.settings.n_jobs = args.cpus

print("Reading data")
adata = ad.read_h5ad(args.input)
df_hvgs = pd.read_pickle(args.hvgs)

print("Integrating data")
if args.method != "unintegrated":
    method = methods[args.method]

    kwargs = {}

    if args.method in ["scvi", "scanvi"]:
        kwargs["save_model"] = True
    if args.method == "scanvi" and args.scvi_model is not None:
        kwargs["scvi_model_path"] = args.scvi_model
    if args.method in ["harmony", "scanorama", "trvaep", "scgen", "scvi", "scanvi", "mnn", "bbknn"]:
        hvg_list = df_hvgs[df_hvgs["highly_variable"]].index.to_list()
        kwargs["hvg"] = hvg_list
    if args.method == "scanvi":
        kwargs["labels_output"] = "scanvi_labels.pkl"
    if args.method == "desc":
        kwargs["ncores"] = args.cpus

        # Check if GPU is available
        if torch.cuda.is_available():
            kwargs["use_gpu"] = True
            # Select random GPU
            kwargs["gpu_id"] = torch.cuda.current_device()

    if args.method in ["scgen", "scanvi"]:
        adata_integrated = method(adata, "batch", "cell_type", **kwargs)
    else:
        adata_integrated = method(adata, "batch", **kwargs)

else:
    adata_integrated = adata.copy()

print("Integration done")

if "X_emb" not in adata_integrated.obsm:
    sc.pp.pca(adata_integrated)
    adata_integrated.obsm["X_emb"] = adata_integrated.obsm["X_pca"]

adata.obsm["X_emb"] = adata_integrated.obsm["X_emb"]

print("Saving integrated data")
adata.write_h5ad(args.output)