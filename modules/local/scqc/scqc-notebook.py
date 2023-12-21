# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
from nxfvars import nxfvars
import matplotlib.pyplot as plt
import seaborn as sns
from qc_plots import plot_qc_metrics, get_stats_df
import pandas as pd

# %%
dataset_id = nxfvars["dataset_id"]
input_adata = nxfvars["input_adata"]
output_adata = nxfvars["output_adata"]
output_stats = nxfvars["output_stats"]
thresholds = {
    key: float(nxfvars[key])
    for key in [
        "min_genes",
        "max_genes",
        "min_counts",
        "max_counts",
        "max_pct_mito"]
}

# %%
adata = sc.read_h5ad(input_adata)

# %%
if adata.__dict__["_raw"] and "_index" in adata.__dict__["_raw"].__dict__["_var"]:
    adata.__dict__["_raw"].__dict__["_var"] = (
        adata.__dict__["_raw"].__dict__["_var"].rename(columns={"_index": "features"})
    )

# %%
# Add fake sample if its not in obs
if "sample" not in adata.obs.columns:
    no_sample = True
    adata.obs["sample"] = "1"
else:
    no_sample = False

# %%
if "mito" not in adata.var.columns:
    adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")

# %%
sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# %%
stats_before = get_stats_df(adata, dataset_id).assign(status=["before_qc"])


# %%
adata.obs.columns

# %%
figwidth = min(max(adata.obs["sample"].unique().size * 0.5, 2), 20)

# %%
fig, ax = plt.subplots(1, 1, figsize=(figwidth, 5))
sc.pl.violin(
    adata, "total_counts", groupby="sample", rotation=90, log=True, cut=0, ax=ax
)

# %%
fig, ax = plt.subplots(1, 1, figsize=(figwidth, 5))
sc.pl.violin(adata, "pct_counts_mito", groupby="sample", rotation=90, ax=ax)

# %%
plot_qc_metrics(adata, **thresholds)

# %%
plot_qc_metrics(adata, cumulative=True, **thresholds)

# %%
# very basic gene filtering - genes with 0 cells cause some downstream processes to fail.
print("Filtering genes")
print(f"    Before: {adata.shape[1]}")
sc.pp.filter_genes(adata, min_counts=3)
print(f"    After: {adata.shape[1]}")

# %%
# Apply thresholds
print("Filter by min_counts")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, min_counts=thresholds["min_counts"])
print(f"    After: {adata.shape[0]}")


print("Filter by max_counts")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, max_counts=thresholds["max_counts"])
print(f"    After: {adata.shape[0]}")


print("Filter by min_genes")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, min_genes=thresholds["min_genes"])
print(f"    After: {adata.shape[0]}")


print("Filter by max_genes")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, max_genes=thresholds["max_genes"])
print(f"    After: {adata.shape[0]}")


print("Filter by max_pct_mito")
print(f"    Before: {adata.shape[0]}")
adata = adata[adata.obs["pct_counts_mito"] < thresholds["max_pct_mito"]].copy()
print(f"    After: {adata.shape[0]}")

# %% [markdown]
# ## After filtering

# %%
plot_qc_metrics(adata, **thresholds)

# %%
plot_qc_metrics(adata, cumulative=True, **thresholds)

# %% [markdown]
# ### Save AnnData object

# %%
if no_sample:
    del adata.obs["sample"]

# %%
stats_after = get_stats_df(adata, dataset_id).assign(status=["after_qc"])
pd.concat([stats_before, stats_after]).to_csv(output_stats, sep="\t", index=False)
adata.write_h5ad(output_adata)

# %%
