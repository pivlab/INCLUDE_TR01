# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Modules

# %%
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
import seaborn as sns
sns.set_theme()

# %%
# where am I?
current_path = Path(os.getcwd()).resolve()
display(current_path)

# %% [markdown]
# # Load data

# %%
pdata = pd.read_pickle(current_path / "plier_model_Z.pkl")

# %%
pdata.shape

# %%
pdata.head()

# %%
dpdata = pd.read_pickle(current_path / "delayed_model_Z.pkl")

# %%
dpdata.shape

# %%
dpdata.head()

# %%
assert pdata.shape == dpdata.shape

# %% [markdown]
# # First attempt: match one-to-one and correlation

# %%
pdata.corrwith(dpdata)


# %% [markdown]
# **Conclusion:** from these correlation numbers, LVs clearly do not match one to one

# %% [markdown]
# # Second attempt: for each LV in model 1, get maximum correlation across all other LVs in model 2

# %%
# Compute the cross-correlation matrix
def cross_correlation(df1, df2, method="pearson"):
    if method == "spearman":
        # compute ranks
        df1 = df1.rank(axis=0)
        df2 = df2.rank(axis=0)
    
    # Normalize each DataFrame (subtract mean and divide by std deviation)
    df1_normalized = (df1 - df1.mean()) / df1.std()
    df2_normalized = (df2 - df2.mean()) / df2.std()

    # Calculate the number of observations (rows)
    n = df1.shape[0]

    # Compute the cross-correlation matrix
    cross_corr_matrix = df1_normalized.T.dot(df2_normalized) / (n - 1)

    return cross_corr_matrix


# %% [markdown]
# ## Using Pearson correlation

# %% [markdown]
# Using the Pearson or Spearman correlation compares two LVs across all genes.

# %%
# make sure function works
_t0 = cross_correlation(pdata, pdata)
_t1 = pdata.corr()
assert np.allclose(_t0, _t1)

# %%
cross_corr_matrix = cross_correlation(pdata, dpdata)

# %%
cross_corr_matrix.shape

# %%
cross_corr_matrix.head()

# %% [markdown]
# ### Get maximum correlation for each LV in PLIER

# %%
plier_lvs_max_corr = cross_corr_matrix.max(axis=1)
display(plier_lvs_max_corr)

# %%
plier_lvs_max_corr.describe()

# %% [markdown]
# **Conclusion:** although half of LVs in PLIER have a correlation larger than 0.89 with LVs in DelayedPLIER, it is clear that these two models are not very similar. If they were identical or very similar, the above numbers should all be close to 1.0:

# %%
cross_correlation(pdata, pdata).max(axis=1).describe()

# %% [markdown]
# Below, I try to identify which LVs are correlated between PLIER and DelayedPLIER:

# %%
plier_lvs_max_corr.sort_values()

# %%
cross_corr_matrix.loc["V233"].sort_values()

# %% [markdown]
# LV 233 in PLIER is highly correlated with LV 363 in DelayedPLIER:

# %%
pdata["V233"].sort_values(ascending=False).head(20)

# %%
dpdata["V363"].sort_values(ascending=False).head(20)

# %% [markdown]
# ### Get maximum correlation for each LV in DelayedPLIER

# %%
dplier_lvs_max_corr = cross_corr_matrix.max(axis=0)
display(dplier_lvs_max_corr)

# %%
dplier_lvs_max_corr.describe()

# %% [markdown]
# ## Using Spearman correlation

# %% [markdown]
# Spearman correlation is more robust to outliers than Pearson.

# %%
# make sure function works
_t0 = cross_correlation(pdata, pdata, method="spearman")
_t1 = pdata.corr(method="spearman")
assert np.allclose(_t0, _t1)

# %%
cross_corr_matrix = cross_correlation(pdata, dpdata, method="spearman")

# %%
cross_corr_matrix.shape

# %%
cross_corr_matrix.head()

# %% [markdown]
# ### Get maximum correlation for each LV in PLIER

# %%
plier_lvs_max_corr = cross_corr_matrix.max(axis=1)
display(plier_lvs_max_corr)

# %%
plier_lvs_max_corr.describe()

# %% [markdown]
# ### Get maximum correlation for each LV in DelayedPLIER

# %%
dplier_lvs_max_corr = cross_corr_matrix.max(axis=0)
display(dplier_lvs_max_corr)

# %%
dplier_lvs_max_corr.describe()

# %% [markdown]
# **Conclusion:** from these, it's even more clear that although some LVs are similar across models, the two models as a whole are not very similar.
# We should expect more similar outputs from PLIER and DelayedPLIER given the same inputs (I know the SVD matrices are not the same, so that's probably the cause of this mismatch between models).

# %% [markdown]
# ## Using Jaccard index

# %% [markdown]
# Here I do the same analysis, but using the Jaccard index and taking the top 1% of genes for each LV.
# The Jaccard index is very simple: it computes the proportion of the top genes that overlap between two LVs.
# This is in a way less noisy than using correlation, because it only takes the top genes of the LVs (which are the most important ones).

# %% [markdown]
# ### Take the top 1% of genes for each LV

# %%
pdata_top_genes = pdata.apply(lambda x: x > x.quantile(0.99))

# %%
pdata_top_genes

# %%
dpdata_top_genes = dpdata.apply(lambda x: x > x.quantile(0.99))

# %%
dpdata_top_genes

# %% [markdown]
# ### Compute cross-Jaccard matrix

# %%
# testing to make sure it's correct
_tmpdf = pd.DataFrame(
    data=1 - pairwise_distances(
        X=pdata_top_genes.iloc[:, :10].T.to_numpy(),
        Y=pdata_top_genes.iloc[:, :10].T.to_numpy(),
        metric="jaccard",
        n_jobs=1
    ),
    index=pdata_top_genes.columns.tolist()[:10],
    columns=dpdata_top_genes.columns.tolist()[:10],
)
_tmp = np.unique(np.diag(_tmpdf))
assert _tmp.shape[0] == 1
assert _tmp[0] == 1.0

assert (_tmpdf.describe().loc["min"] >= 0.0).all()
assert (_tmpdf.describe().loc["max"] <= 1.0).all()

# %%
cross_corr_matrix = pd.DataFrame(
    data=1 - pairwise_distances(
        X=pdata_top_genes.T.to_numpy(),
        Y=dpdata_top_genes.T.to_numpy(),
        metric="jaccard",
        n_jobs=1
    ),
    index=pdata_top_genes.columns.tolist(),
    columns=dpdata_top_genes.columns.tolist(),
)

# %%
cross_corr_matrix.shape

# %%
cross_corr_matrix.head()

# %% [markdown]
# ### Get maximum Jaccard for each LV in PLIER

# %%
plier_lvs_max_corr = cross_corr_matrix.max(axis=1)
display(plier_lvs_max_corr)

# %%
plier_lvs_max_corr.describe()

# %%
plier_lvs_max_corr.sort_values().tail(10)

# %% [markdown]
# **Conclusion:** the above list does not show LV 233 from PLIER, which was very correlated with LV 363 from DelayedPLIER when using Pearson.
# In the above list, we see other LVs more "similar" (according to Jaccard), because they share more genes at the top 1%.
#
# Let's see what happens with LV 233 in PLIER and LV 363 in DelayedPLIER:

# %%
cross_corr_matrix.loc["V233"].sort_values()

# %%
pdata["V233"].sort_values().tail(67).index.intersection(
    dpdata["V363"].sort_values().tail(67).index
).shape

# %% [markdown]
# LV 233 and LV 363 share 60 genes out of 67. But LV 78, for example, share all the top 1% genes with LV 78 in DelayedPLIER (so Jaccard is 1.0):

# %%
cross_corr_matrix.loc["V78"].sort_values()

# %%
pdata["V78"].sort_values().tail(67).index.intersection(
    dpdata["V78"].sort_values().tail(67).index
).shape

# %% [markdown]
# ### Get maximum Jaccard for each LV in DelayedPLIER

# %%
dplier_lvs_max_corr = cross_corr_matrix.max(axis=0)
display(dplier_lvs_max_corr)

# %%
dplier_lvs_max_corr.describe()

# %% [markdown]
# # How does it look the correlation between LVs in the PLIER model?

# %% [markdown]
# This part is just to take a look at the correlation matrix. I'm not comparing the two models here (it's clear they differ).

# %%
pdata_corr = pdata.corr(method="spearman")

# %%
pdata_corr.shape

# %%
pdata_corr

# %%
# taken from here: https://seaborn.pydata.org/examples/structured_heatmap.html
g = sns.clustermap(
    pdata_corr,
    center=0,
    cmap="vlag",
    # row_colors=network_colors,
    # col_colors=network_colors,
    dendrogram_ratio=(.1, .2),
    cbar_pos=(.02, .32, .03, .2),
    # linewidths=.75, figsize=(12, 13)
    yticklabels=False,
    xticklabels=False,
)

g.ax_row_dendrogram.remove()

# %%
