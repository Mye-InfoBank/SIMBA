# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars

sc.settings.set_figure_params(figsize=(5, 5))
from pathlib import Path
from scanpy_helpers.annotation import AnnotationHelper
import progeny
import dorothea
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
import altair as alt
import re
import statsmodels.stats
import warnings
import scipy.stats
from tqdm.auto import tqdm
import numpy as np
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import itertools
from operator import or_
from functools import reduce
import scanpy_helpers as sh

# %%
threadpool_limits(32)

# %%
ah = AnnotationHelper()

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
# For the patient stratification, treat the two batches of the UKIM-V dataset as one
adata.obs["dataset"] = adata.obs["dataset"].str.replace("UKIM-V-2", "UKIM-V")

# %%
sc.pl.umap(adata, color=["cell_type_major", "origin"], ncols=1)

# %%
adata.obs["cell_type_major"].value_counts().sort_index()

# %%
# # TODO cDC1 is missing; maybe use the new "cell_type_major" annotation
cell_types = {
    "tumor": ["Tumor cells"],
    "healthy epithelial": [
        "Alveolar cell type 2",
        "Alveolar cell type 1",
        "Ciliated",
        "Club",
        "Goblet",
    ],
    "immune": [
        "B cell",
        "DC mature",
        "Macrophage",
        "Macrophage FABP4+",
        "Mast cell",
        "Monocyte",
        "NK cell",
        # "Neutrophils",
        "Plasma cell",
        "T cell CD4",
        "T cell CD8",
        "T cell regulatory",
        "cDC1",
        "cDC2",
        "pDC",
    ],
    "structural": ["Endothelial cell", "Stromal"],
}
assert reduce(or_, [set(v) for v in cell_types.values()]) | {
    "other",
    "Neutrophils",
} == set(adata.obs["cell_type_major"]), "Some cell-types are missing"

# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
    :,
]

# %%
adata_tumor_cells = adata_primary_tumor[
    adata_primary_tumor.obs["cell_type_major"] == "Tumor cells",
    :,
]

# %% [markdown]
# ## Tumor subtypes

# %%
ad_tumor_subtypes = sc.AnnData(
    X=adata_tumor_cells.obs.groupby(["patient", "cell_type_tumor"], observed=True)
    .size()
    .reset_index(name="n")
    .assign(
        cell_type_tumor=lambda x: x["cell_type_tumor"]
        .str.replace("Tumor cells ", "")
        .str.replace(" mitotic", "")
    )
    .pivot_table(values="n", columns="cell_type_tumor", index="patient", fill_value=0),
)
ad_tumor_subtypes.obs = ad_tumor_subtypes.obs.join(
    adata_tumor_cells.obs.loc[:, ["patient", "condition"]]
    .drop_duplicates()
    .set_index("patient")
)


# %%
sc.pp.normalize_total(ad_tumor_subtypes, target_sum=1)

# %%
sc.pl.heatmap(
    ad_tumor_subtypes,
    groupby="condition",
    var_names=ad_tumor_subtypes.var_names,
    swap_axes=True,
    figsize=(14, 1.5),
)

# %%
ad_tumor_subtypes.obs["predominant_tumor_subtype"] = ad_tumor_subtypes.var_names[
    np.argmax(ad_tumor_subtypes.X, axis=1)
]

# %%
ad_tumor_subtypes.obs


# %% [markdown]
# ## Infiltration patterns

# %%
def get_cell_type_group(ct):
    for group, cts in cell_types.items():
        if ct in cts:
            return group
    return "other"


# %%
major_cell_types_df = (
    (
        adata_primary_tumor.obs.assign(
            cell_type_group=lambda x: x["cell_type_major"].apply(get_cell_type_group)
        )
        .groupby(["dataset", "patient", "cell_type_group"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type_group",
            index=["dataset", "patient"],
            fill_value=0,
        )
        .drop("other", axis="columns")
    )
    .assign(immune_tumor_ratio=lambda x: np.log2((x["immune"] + 1) / (x["tumor"] + 1)))
    .assign(
        structural_tumor_ratio=lambda x: np.log2(
            (x["structural"] + 1) / (x["tumor"] + 1)
        )
    )
)


# %%
ad_ti_ratio = sc.AnnData(
    obs=major_cell_types_df.drop(
        columns=["immune_tumor_ratio", "structural_tumor_ratio"]
    ),
    X=major_cell_types_df.loc[:, ["immune_tumor_ratio", "structural_tumor_ratio"]],
)

# %%
ad_ti_ratio.obs = ad_ti_ratio.obs.reset_index().set_index("patient")

# %%
sc.pl.matrixplot(
    ad_ti_ratio,
    groupby="patient",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7
    # vcenter=0
)

# %%
sc.pp.regress_out(ad_ti_ratio, "dataset")

# %%
ad_ti_ratio.obs.index.name = "index"
ad_ti_ratio.obs["patient"] = ad_ti_ratio.obs.index

# %%
sc.pl.matrixplot(
    ad_ti_ratio,
    groupby="patient",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7,
    # vcenter=0
)

# %% [markdown]
# ## All cell-types

# %%
ad_cts = sc.AnnData(
    X=(
        adata_primary_tumor.obs.assign(
            cell_type_group=lambda x: x["cell_type_major"].apply(get_cell_type_group)
        )
        .loc[lambda x: x["cell_type_group"] != "other", :]
        .groupby(["dataset", "patient", "cell_type"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type",
            index=["dataset", "patient"],
            fill_value=0,
        )
    )
)
ad_cts.obs = ad_cts.obs.reset_index().set_index("patient")

# %%
sc.pp.normalize_total(ad_cts, target_sum=1)

# %%
sc.pl.matrixplot(
    ad_cts,
    var_names=ad_cts.var_names,
    groupby="patient",
    dendrogram=False,
    swap_axes=True,
    cmap="viridis",
    # vmin=-0.25,
    # vmax=0.25
    # # vmin=0,
    # vmax=1,
    standard_scale="var",
)

# %%
sc.pp.regress_out(ad_cts, "dataset")

# %%
ad_cts.obs["patient"] = ad_cts.obs.index
ad_cts.obs.index.name = "index"

# %%
sc.tl.dendrogram(ad_cts, groupby="patient", use_rep="X", optimal_ordering=True)

# %%
sc.pl.matrixplot(
    ad_cts,
    var_names=ad_cts.var_names,
    groupby="patient",
    dendrogram=True,
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %% [markdown]
# ## immune cell patterns
# (only immune cells)

# %%
ad_immune = sc.AnnData(
    X=(
        adata_primary_tumor.obs.assign(
            cell_type_group=lambda x: x["cell_type_major"].apply(get_cell_type_group)
        )
        .loc[lambda x: x["cell_type_group"] == "immune", :]
        .groupby(["dataset", "patient", "cell_type"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type",
            index=["dataset", "patient"],
            fill_value=0,
        )
    )
)
ad_immune.obs = ad_immune.obs.reset_index().set_index("patient")

# %%
sc.pp.normalize_total(ad_immune, target_sum=1)

# %%
sc.pl.matrixplot(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="patient",
    dendrogram=False,
    swap_axes=True,
    cmap="viridis",
    # vmin=-0.25,
    # vmax=0.25
    # # vmin=0,
    # vmax=1,
    standard_scale="var",
)

# %%
sc.pp.regress_out(ad_immune, "dataset")

# %%
ad_immune.obs["patient"] = ad_immune.obs.index
ad_immune.obs.index.name = "index"

# %%
sc.tl.dendrogram(ad_immune, groupby="patient", use_rep="X", optimal_ordering=True)

# %%
sc.pl.matrixplot(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="patient",
    dendrogram=True,
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %% [markdown]
# ## Clustering and stratification

# %% [markdown]
# ### Immune/structural ratio

# %%
imm = ad_ti_ratio[:, "immune_tumor_ratio"].X > 0
structural = ad_ti_ratio[:, "structural_tumor_ratio"].X > 0


def get_state(i, s):
    res = []
    res.append("I" if i else "-")
    res.append("S" if s else "-")
    return "/".join(res)


ad_ti_ratio.obs["group"] = [get_state(i, s) for i, s in zip(imm, structural)]

# %%
sc.pl.heatmap(
    ad_ti_ratio,
    groupby="group",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7
    # vcenter=0
)

# %% [markdown]
# ### immune-cell infiltration patterns

# %%
ad_immune.obs_names = ad_immune.obs_names.values.astype(str)

# %%
ad_immune_sub = ad_immune[
    list(
        set(
            ad_ti_ratio[
                ad_ti_ratio.obs["group"].str.contains("I", regex=False), :
            ].obs_names
        )
        & set(ad_immune.obs_names)
    ),
    :,
]

# %%
sc.pp.neighbors(ad_immune_sub, use_rep="X")

# %%
sc.tl.leiden(ad_immune_sub, resolution=0.75)

# %%
sc.pl.heatmap(
    ad_immune_sub,
    var_names=ad_immune.var_names,
    groupby="leiden",
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
ad_immune_sub.obs["immune_type"] = [
    {"0": "T", "1": "B", "2": "M"}.get(x, "") for x in ad_immune_sub.obs["leiden"]
]

# %%
sc.pl.heatmap(
    ad_immune_sub,
    var_names=ad_immune.var_names,
    groupby="immune_type",
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %% [markdown]
# ## Make figure

# %%
# Patients with "primary tumor" samples
adata_primary_tumor.obs["patient"].nunique()

# %%
# Minus 4 patients with no "tumor cells" in primary tumor samplese
adata_tumor_cells.obs["patient"].nunique()

# %%
ad_tumor_subtypes.obs

# %%
ad_ti_ratio.obs["patient"].nunique()


# %%
def get_stratum(infiltration_group, immune_type):
    stratum = ["-", "-"]
    if "I" in infiltration_group:
        stratum[0] = immune_type
    if "S" in infiltration_group:
        stratum[1] = "S"
    return "/".join(stratum)


# %%
adata.obs.columns

# %%
plot_df = (
    adata.obs.loc[
        :,
        [
            "patient",
            "sex",
            "uicc_stage",
            "ever_smoker",
            "age",
            "tumor_stage",
            "dataset",
            "platform",
        ],
    ]
    .drop_duplicates()
    .set_index("patient")
    .join(ad_tumor_subtypes.obs, how="inner")
    .join(ad_ti_ratio.obs.loc[:, ["group"]])
)
plot_df["immune_type"] = ad_immune_sub.obs["immune_type"].astype(str)
plot_df.loc[plot_df["immune_type"].isnull(), "immune_type"] = "-"
plot_df["stratum"] = [
    get_stratum(infil, imm)
    for infil, imm in zip(plot_df["group"], plot_df["immune_type"])
]
remap_tumor_type = {
    "PPC": "NSCLC",
    "LUSC": "LSCC",
    "LCLC": "NSCLC",
}
plot_df["condition"] = [remap_tumor_type.get(x, x) for x in plot_df["condition"]]
plot_df["sex"] = plot_df["sex"].cat.add_categories("unknown").fillna("unknown")
plot_df.rename(
    columns={
        "condition": "tumor_type_annotated",
        "predominant_tumor_subtype": "tumor_type_inferred",
        "group": "infiltration_state",
        "immune_type": "immune_infiltration",
        "stratum": "TMIG",  # tumor micro environment stratum
    },
    inplace=True,
)
plot_df.sort_values(["TMIG", "tumor_type_inferred"], inplace=True)

# %%
plot_df = plot_df.reset_index()

# %%
plot_df


# %%
def get_row(col, color_scale=None):
    if color_scale is None:
        color_scale = col
    return (
        alt.Chart(plot_df.assign(ylab=col))
        .mark_rect()
        .encode(
            x=alt.X(
                "patient:N",
                axis=alt.Axis(labels=False, ticks=False, title=None),
                sort=plot_df["patient"].values,
            ),
            y=alt.Y("ylab", axis=alt.Axis(title=None)),
            color=alt.Color(
                col,
                scale=sh.colors.altair_scale(color_scale),
                legend=alt.Legend(columns=4),
            ),
        )
        .properties(width=800)
    )


p0 = (
    alt.vconcat(
        get_row("tumor_type_annotated", "condition")
        & get_row("tumor_type_inferred", "condition"),
        get_row("sex"),
        get_row("tumor_stage"),
        get_row("dataset"),
        get_row("platform"),
        # get_row("infiltration_state"),
        # get_row("immune_infiltration"),
        get_row("immune_infiltration"),
    ).resolve_scale(color="independent")
    # .resolve_legend("shared")
    .configure_concat(spacing=0)
)
p0


# %%
def scale_range(a):
    return a / max(np.abs(np.max(a)), np.abs(np.min(a)))


def scale_01(a):
    return (a - np.min(a)) / np.ptp(a)


# %%
tmp_ad = ad_ti_ratio[plot_df["patient"], ["immune_tumor_ratio"]]
heatmap_df = (
    pd.DataFrame(
        scale_range(tmp_ad.X), columns=tmp_ad.var_names, index=tmp_ad.obs_names
    )
    .reset_index()
    .rename(columns={"index": "patient"})
    .melt(id_vars="patient")
)
p1 = (
    alt.Chart(heatmap_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N",
            sort=plot_df["patient"].values,
            axis=alt.Axis(ticks=False, labels=False, title=None),
        ),
        y=alt.Y("cell_type_group:N", axis=alt.Axis(title=None)),
        color=alt.Color(
            "value", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-0.9, 0.9])
        ),
    )
    .properties(width=800, height=20)
)

# %%
tmp_ad = ad_immune[
    plot_df.loc[plot_df["patient"].isin(ad_immune.obs_names), "patient"], :
].copy()
# scale by tumor immune ratio
tmp_ad.X = scale_range(tmp_ad.X) * scale_01(
    ad_ti_ratio[plot_df["patient"], "immune_tumor_ratio"].X
)

heatmap_df = (
    pd.DataFrame(tmp_ad.X, columns=tmp_ad.var_names, index=tmp_ad.obs_names)
    .reset_index()
    .rename(columns={"index": "patient"})
    .melt(id_vars="patient")
)
p2 = (
    alt.Chart(heatmap_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N",
            sort=plot_df["patient"].values,
            axis=alt.Axis(ticks=False, labels=False),
        ),
        y=alt.Y(
            "cell_type:N",
            sort=[
                "B cell",
                "Plasma cell",
                "Mast cell",
                "Macrophage",
                "Macrophage FABP4+",
                "Monocyte",
                "NK cell",
                "T cell CD4",
                "T cell CD8",
                "T cell regulatory",
                "DC mature",
                "cDC2",
                "pDC",
            ],
            axis=alt.Axis(title=None),
        ),
        color=alt.Color("value", scale=alt.Scale(scheme="redblue")),
    )
    .properties(width=800, height=120)
)

# %%
p0 & (p1 & p2).resolve_scale(x="shared")

# %%
tmp_ad.shape

# %% [markdown]
# ## groups by histological subtype

# %%
tumor_type_total = (
    plot_df.groupby(["tumor_type_annotated"]).size().reset_index(name="total")
)

# %%
tmp_df = (
    plot_df.groupby(["immune_infiltration", "tumor_type_annotated"])
    .size()
    .reset_index(name="n")
    .merge(tumor_type_total)
    .assign(frac_of_total=lambda x: x["n"] / x["total"])
)

# %%
alt.Chart(tmp_df).mark_bar().encode(
    x=alt.X("tumor_type_annotated", title=None),
    y="frac_of_total",
    color=alt.Color("tumor_type_annotated", scale=sh.colors.altair_scale("tumor_type")),
).facet(column="immune_infiltration")

# %% [markdown]
# # Write output

# %%
plot_df.to_csv(
    "{}/patient_stratification.csv".format(
        nxfvars.get("artifact_dir", "/home/sturm/Downloads")
    )
)

# %%
