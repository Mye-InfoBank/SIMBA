# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

from nxfvars import nxfvars
import scanpy as sc
import numpy as np
import itertools
from tqdm import tqdm
import os

# %%
paths = [file for file in os.listdir(".") if file.endswith(".h5ad")]

if len(paths) != 1:
    raise ValueError("Expected only one file in this directory, but found: {}".format(paths))
adata_path = paths[0]
adata_path


# %%
adata = sc.read_h5ad(adata_path)
adata

# %%
sc.pp.pca(adata)

# %%
sc.external.pp.harmony_integrate(adata, key='batch')

# %%
os.makedirs("artifacts", exist_ok=True)

adata.write_h5ad("artifacts/harmony.h5ad")