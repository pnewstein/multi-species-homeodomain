"""
makes dotplots for bipolar cells and lamina cells reading h5ad files
"""

from pathlib import Path
import sys
from typing import Sequence

import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import seaborn as sns


try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()


IMG_DIR = HERE / "../imgs"


def get_expressed_htdf(
    adata, transcript_thresh: float, hox_genes: pd.DataFrame, varname: str
) -> list[str]:
    hdtf_dict = dict(
        zip(hox_genes["external_gene_name"].values, hox_genes["is_hdtf"].values)
    )
    out: set[str] = set()
    missing_genes: list[str] = []
    for cluster_name, sub_df in adata.obs.groupby(varname):
        sub_adata = adata[sub_df.index, :]
        mean_expression_array_like = sub_adata[sub_df.index, :].X.mean(axis=0)
        if isinstance(mean_expression_array_like, np.matrix):
            mean_expression_array = np.array(mean_expression_array_like).flatten()
        else:
            mean_expression_array = mean_expression_array_like
        is_expressed = mean_expression_array > transcript_thresh
        expressed_genes = adata.var_names[is_expressed]
        # drop mitochondrial
        expressed_genes = expressed_genes[~expressed_genes.str.startswith("mt-")]
        # drop others
        expressed_genes = expressed_genes.drop(
            [
                "E130218I03Rik",
                "Fam19a3",
                "Gucy1b3",
                "Malat1",
                "Meg3",
                "Mir124-2hg",
                "Pnmal2",
                "Xist",
                "Atpif1",
                "BC030499",
                "Hist3h2a",
                "Skp1a",
                "2010107E04Rik",
                "Atp5b",
                "Sept7",
                "Atp5f1",
                "Atp5g3",
                "Gucy1a3",
                "A730046J19Rik",
                "Gm4792",
                "Gas5",
                "Gnb2l1",
                "Lect1",
                "H2afy",
                "Atp5c1",
                "Atp5j",
                "Gm37583",
                "Hist3h2ba",
            ],
            errors="ignore",
        )
        found_bools = np.isin(expressed_genes, hox_genes["external_gene_name"])
        missing_genes.extend(list(expressed_genes[~found_bools]))
        expressed_genes = expressed_genes[found_bools]
        expressed_hdtfs = expressed_genes[[hdtf_dict[g] for g in expressed_genes]]
        out = out.union(np.array(expressed_hdtfs))
    Path("missing_genes.txt").write_text("\n".join(set(missing_genes)))
    return list(out)


def plot_dotplot(
    adata: sc.AnnData,
    obs_name: str,
    genes: list[str],
    thresh: float,
    cmap: str,
    order: Sequence[str] | None = None,
):
    cluster_names = np.unique(adata.obs[obs_name])
    expression_df = adata[:, np.array(genes)].to_df()
    expression_df[obs_name] = adata.obs[obs_name]
    clusterwise_expression = pd.DataFrame(np.nan, index=cluster_names, columns=genes)
    for cluster, sub_df in expression_df.groupby(obs_name, observed=True):
        if cluster == "other":
            continue
        clusterwise_expression.loc[cluster, :] = sub_df.loc[:, genes].mean()
    if order is None:
        hdtfs_in_order = (clusterwise_expression > thresh).sum().sort_values().index
    else:
        hdtfs_in_order = list(order)
    clusters_in_order = np.sort(clusterwise_expression.index).tolist()
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
    return fig, axs


def bipolar_dots():
    adata = sc.read_h5ad(HERE / "../data/bpc.h5ad")
    # adata.obs["category"] = adata.obs["category"].astype(pd.Categorical(np.unique(adata.obs["category"])).dtype)
    hox_genes = pd.read_csv(HERE / "../data/mouse_hdtf.csv", index_col=0)
    thresh = 0.355
    high_expression_hdtf = get_expressed_htdf(adata, thresh, hox_genes, "cluster")
    plt.style.use(HERE / "paper.mplstyle")
    fig, ax = plot_dotplot(
        adata, "cluster", high_expression_hdtf, thresh, order=None, cmap="plasma_r"
    )
    fig.set_size_inches((8.21, 5.76))
    fig.savefig(IMG_DIR / "bpc_dots.svg")
    fig: Figure
    fig, ax = plt.subplots(constrained_layout=True)  # type: ignore
    sc.pl.umap(adata, color="cluster", ax=ax, legend_loc="on data")
    fig.set_size_inches((2.2, 1.7))
    fig.savefig(IMG_DIR / "bpc_umap.svg")


def main():
    if sys.argv[1] == "bpc":
        bipolar_dots()
    elif sys.argv[1] == "lamina":
        lamina_dots()
    elif sys.argv[1] == "retina":
        retina_dots()
    else:
        raise ValueError("unknown command")


def lamina_dots():
    adata = sc.read_h5ad(HERE / "../data/chundi_lamina.h5ad")
    adata.obs["category"] = adata.obs["category"].astype(pd.Categorical(np.unique(adata.obs["category"])).dtype)
    # first plot a umap
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, spread=3)
    plt.style.use(HERE / "paper.mplstyle")
    sc.pl.umap(adata, color="category")
    fig = plt.gcf()
    fig.set_size_inches((8.45, 7.84))
    fig.savefig("lamina_umap.svg")
    hox_genes = pd.read_csv(HERE / "../data/drosophila_hdtf.csv", index_col=0)
    thresh = .5
    high_expression_hdtf = get_expressed_htdf(adata, thresh, hox_genes, "category")
    order = "ap", "pdm3", "bsh", "zfh1", "zfh2", "onecut", "scro"
    assert set(order) == set(high_expression_hdtf)
    fig, ax = plot_dotplot(adata, "category", high_expression_hdtf, thresh, order=order, cmap="plasma_r")
    fig.set_size_inches((5, 3.5))
    fig.savefig(HERE / "../imgs/lamina_dots.svg")

if __name__ == "__main__":
    main()
