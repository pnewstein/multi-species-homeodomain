"""
makes dotplots for bipolar cells and lamina cells reading h5ad files
"""

from pathlib import Path
import sys
from typing import Sequence
from itertools import chain

import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import seaborn as sns
import umap


try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()

sys.path.append(str(HERE))
from query_domains import get_expressed_hdtf


IMG_DIR = HERE / "../imgs"


def plot_dotplot(
    adata: sc.AnnData,
    obs_name: str,
    thresh: float,
    cmap: str,
    order: Sequence[str] | None = None,
    use_umap=False,
    fly=False,
) -> Figure:
    cluster_hdtf_expression = get_expressed_hdtf(adata, thresh, obs_name, fly)
    genes = list(set(chain.from_iterable(cluster_hdtf_expression.values())))
    cluster_names = np.unique(adata.obs[obs_name])
    expression_df = adata[:, np.array(genes)].to_df()
    expression_df[obs_name] = adata.obs[obs_name]
    clusterwise_expression = pd.DataFrame(np.nan, index=cluster_names, columns=genes)
    for cluster, sub_df in expression_df.groupby(obs_name, observed=True):
        if cluster == "other":
            continue
        clusterwise_expression.loc[cluster, :] = sub_df.loc[:, genes].mean()
    if order is None:
        number_above_thresh = (clusterwise_expression > thresh).sum()
        sum_expression = clusterwise_expression.sum(axis=0)
        hdtfs_in_order = (
            pd.DataFrame([number_above_thresh, sum_expression])
            .T.sort_values([0, 1])
            .index
        )
    else:
        hdtfs_in_order = list(order)
        assert set(genes) == set(hdtfs_in_order), "Wrong genes are plotted"
    clusters_in_order = np.sort(clusterwise_expression.index).tolist()
    if use_umap:
        reduced = (
            umap.UMAP(n_components=1, random_state=100)
            .fit_transform(clusterwise_expression)
            .flatten()
        )
        clusters_in_order = pd.Series(reduced, index=cluster_names).sort_values().index
    ax: Axes
    fig: Figure
    fig, ax = plt.subplots(constrained_layout=True)  # type: ignore
    edge_color = ax.spines["left"].get_edgecolor()
    dp = sc.pl.DotPlot(
        adata,
        var_names=hdtfs_in_order,
        groupby=obs_name,
        categories_order=clusters_in_order,
        ax=ax,
    )
    dp.style(cmap=cmap, dot_edge_color=edge_color)
    dp.swap_axes()
    axs = dp.show(return_axes=True)
    sns.despine(fig)
    return fig


def bipolar_dots():
    adata = sc.read_h5ad(HERE / "../data/bpc.h5ad")
    thresh = 0.355
    plt.style.use(HERE / "paper.mplstyle")
    fig = plot_dotplot(adata, "cluster", thresh, order=None, cmap="plasma_r")
    fig.set_size_inches((8.21, 5.76))
    fig.savefig(IMG_DIR / "bpc_dots.svg")


def retina_dots():
    adata = sc.read_h5ad(HERE / "../data/MRCA_full_normcounts.h5ad")
    adata.obs["majorclass"] = adata.obs["majorclass"].astype(str)
    adata.obs.loc[adata.obs["majorclass"].isin(["Rod", "Cone"]), "majorclass"] = "PR"
    adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")
    size_dict = {
        "AC": (16.61, 7.49),
        "RGC": (12.5, 7.77),
        "BC": (12.85, 7.65),
        "PR": (2.02, 3.83),
        "HC": (1.86, 3.2),
    }
    figures: list[Figure] = []
    for cell_class, size in size_dict.items():
        adata_in_class = adata[adata.obs["majorclass"] == cell_class]
        use_umap = False if cell_class == "PR" else True
        with plt.style.context("paper.mplstyle"):
            fig = plot_dotplot(
                adata_in_class,
                "celltype",
                thresh=0.35,
                cmap="plasma_r",
                use_umap=use_umap,
            )
        fig.set_size_inches(size)
        fig.savefig(HERE / f"imgs/{cell_class}.svg")
        figures.append(fig)


def lamina_dots():
    adata = sc.read_h5ad(HERE / "../data/chundi_lamina.h5ad")
    adata.obs["category"] = adata.obs["category"].astype(
        pd.Categorical(np.unique(adata.obs["category"])).dtype
    )
    # first plot a umap
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, spread=3)
    plt.style.use(HERE / "paper.mplstyle")
    sc.pl.umap(adata, color="category")
    fig = plt.gcf()
    fig.set_size_inches((8.45, 7.84))
    fig.savefig("lamina_umap.svg")
    hox_genes = pd.read_csv(HERE / "../data/drosophila_hdtf.csv", index_col=0)
    thresh = 0.5
    high_expression_hdtf = get_expressed_htdf(adata, thresh, hox_genes, "category")
    order = "ap", "pdm3", "bsh", "zfh1", "zfh2", "onecut", "scro"
    assert set(order) == set(high_expression_hdtf)
    fig = plot_dotplot(
        adata, "category", thresh, order=order, cmap="plasma_r"
    )
    fig.set_size_inches((5, 3.5))
    fig.savefig(HERE / "../imgs/lamina_dots.svg")


def main():
    if sys.argv[1] == "bpc":
        bipolar_dots()
    elif sys.argv[1] == "lamina":
        lamina_dots()
    elif sys.argv[1] == "retina":
        retina_dots()
    else:
        raise ValueError("unknown command")


if __name__ == "__main__":
    main()
