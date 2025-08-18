"""
Load quantified insitu data and compare to scRNAseq
"""

from pathlib import Path

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import scanpy as sc

COLORS = [
    (0.0, 0.5, 1.0),
    (0.9387206358364905, 0.5058556114436428, 0.0008967132642730968),
    (0.9969681618399396, 0.0025991822706286083, 0.735197084889194),
    (0.02117113131977666, 0.9237412617513718, 0.5062751061006804),
    (0.06509019506991442, 0.04143205627466484, 0.76589626060724),
    (0.04356896683344247, 0.4467629465350331, 0.10443121077210948),
    (0.5257688103789879, 0.6481491483812926, 0.9920661809390425),
]

try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()

COL_PREFIX = "eroded"
# COL_PREFIX = "full"

BPC_MARKER_THRESHS = {"Vsx2": 1, "Otx2": 3}
bpc_marker_threshs = {"Vsx2":1e-4, "Otx2":4e-5}

gene_threshs = {
    "Meis2": 0.2,
    "Lhx3": 0.2,
    "Vsx2": 0.4,
    "Zfhx4": 0.17,
    "Lhx4": 0.3,
    "Sebox": 0.27,
    "Otx2": 0.04,
    "Zeb2": 0.17,
    "Vsx1": 0.3,
    "Six3": 0.4,
    "Irx5": 0.3,
    "Isl1": 0.15,
}


def get_in_situ_csvs(
    paths_chan_dict: dict[Path, tuple[str, str, str]] | None = None
) -> list[pd.DataFrame]:
    if paths_chan_dict is None:
        paths_chan_dict = {
            HERE
            / "../quantified_image_data/Well1_405WGA_488Lhx4_546Irx5_647Otx2.csv": (
                "Lhx4",
                "Irx5",
                "Otx2",
            ),
            HERE
            / "../quantified_image_data/Well2_405WGA_488Vsx2_546Vsx1_647Lhx4.csv": (
                "Vsx2",
                "Vsx1",
                "Lhx4",
            ),
            HERE
            / "../quantified_image_data/Well2_405wga_488Vsx2_546Lhx3_647Six3-1.csv": (
                "Vsx2",
                "Lhx3",
                "Six3",
            ),
            HERE
            / "../quantified_image_data/Well3_405WGA_488Zeb2_546Lhx3_647Sebox.csv": (
                "Zeb2",
                "Lhx3",
                "Sebox",
            ),
            HERE
            / "../quantified_image_data/Well3_405wga_488Zeb2_546Otx2_647Zfhx4-1.csv": (
                "Zeb2",
                "Otx2",
                "Zfhx4",
            ),
            HERE
            / "../quantified_image_data/Well4_405WGA_488Sebox_546Vsx2_647Isl1.csv": (
                "Sebox",
                "Vsx2",
                "Isl1",
            ),
        }
    out: list[pd.DataFrame] = []
    for path, chan_names in paths_chan_dict.items():
        in_df = pd.read_csv(path, index_col=0)
        out_series = [
            pd.Series(in_df[f"{COL_PREFIX}_size"], name="Cell size", index=in_df.index)
        ]
        for i, gene in enumerate(chan_names, 1):
            out_series.append(
                pd.Series(in_df[f"{COL_PREFIX}-{i}"], name=gene, index=in_df.index)
            )
        out_df = pd.DataFrame(out_series).T.drop(0)
        out.append(out_df)
    return out


def preproc(in_df: pd.DataFrame, filter_bpc=True) -> pd.DataFrame:
    min_pix = 10_000
    max_pix = 50_000
    # filter out wrong size
    df = in_df.loc[(min_pix < in_df["Cell size"]) & (in_df["Cell size"] < max_pix)]
    # filter out non bpc
    gene_columns = in_df.columns[1:]
    try:
        bpc_marker = next(c for c in gene_columns if c in BPC_MARKER_THRESHS)
    except StopIteration:
        raise ValueError("No bpc marker")
    if filter_bpc:
        bpc_df = df.loc[df[bpc_marker] > BPC_MARKER_THRESHS[bpc_marker]]
    else:
        bpc_df = df
    # log scale norm by cell size
    out_df = pd.DataFrame(np.nan, columns=gene_columns, index=bpc_df.index)
    for col in gene_columns:
        out_df[col] = np.log1p(bpc_df[col]) / bpc_df["Cell size"].astype(float)
    out_df["Cell size"] = bpc_df["Cell size"]
    return out_df


def plot_in_situ_on_axs(
    in_situ_df: pd.DataFrame,
    sc_df: pd.DataFrame,
    ax: Axes,
    scale: float,
    debug_plot=False,
):
    in_situ_df = in_situ_df.copy()
    if set(in_situ_df.columns) == set(["Vsx2", "Meis2", "Otx2", "Cell size"]):
        exp_genes = ["Meis2", "Vsx2"]
    else:
        exp_genes = [
            c
            for c in in_situ_df.columns
            if c not in list(BPC_MARKER_THRESHS.keys()) + ["Cell size"]
        ]
    gene_x, gene_y = exp_genes
    in_situ_df.loc[:, exp_genes] = scale * in_situ_df.loc[:, exp_genes]
    sns.scatterplot(
        sc_df, x=gene_x, y=gene_y, ax=ax, hue=tuple(exp_genes), alpha=1, s=10
    )
    for handle in ax.legend().legend_handles:
        assert handle is not None
        handle.set_alpha(1)
    sns.scatterplot(in_situ_df, x=gene_x, y=gene_y, ax=ax, alpha=0.5, c="k")
    clusters = np.unique(sc_df["cluster"])
    assert len(clusters) < 15
    if debug_plot:
        _, axs = plt.subplots(3, 5, sharex=True, sharey=True, constrained_layout=True)
        assert isinstance(axs, np.ndarray)
        for ax, (cluster, df) in zip(
            axs.ravel(), sc_df.groupby("cluster", observed=True)
        ):
            assert isinstance(ax, Axes)
            ax.set_title(str(cluster))
            sns.scatterplot(df, x=gene_x, y=gene_y, ax=ax, alpha=0.5)


combo_regex = {
    "BC5-8": r"BC[5678)].*",
    "RBC": "RBC",
    "BC1": "BC1[AB]",
    "BC2/6": "BC[26]",
    "BC8/9": "BC8/9",
    "BC1B": "BC1B",
    "BC2+BC3A": r"BC[23]A?",
    "BC4/5+BC3B": r"BC[45].?|BC3B",
    "BC5D": "BC5D",
    "BC3B+BC5A": "BC3B|BC5A",
    "BC7": "BC7",
}

combo_dicts = {
    ("Sebox", "Isl1"): ("BC5-8", "RBC"),
    ("Zeb2", "Zfhx4"): ("BC1", "BC8/9"),
    ("Lhx3", "Six3"): ("BC2/6", "BC1B"),
    ("Vsx1", "Lhx4"): ("BC2+BC3A", "BC7", "BC4/5+BC3B"),
    ("Meis2", "Vsx2"): ("BC5D", "BC3B+BC5A"),
    ("Lhx4", "Irx5"): ("BC2+BC3A", "BC4/5+BC3B"),
}


def average_sc_rnaseq(df: pd.DataFrame, n_examples: int, n_sample: int) -> pd.DataFrame:
    """
    Make n_examples different estimentes of the average
    """
    out_series: list[pd.Series] = []
    for cluster, sub_df in df.groupby("cluster", observed=False):
        del sub_df["cluster"]
        for _ in range(n_examples):
            ser = sub_df.sample(n_sample).mean()
            assert isinstance(ser, pd.Series)
            ser["cluster"] = cluster
            out_series.append(ser)
    return pd.DataFrame(out_series)


def combine_groups(df: pd.DataFrame) -> pd.DataFrame:
    "combines groups together according to combo_dicts"
    df = df.copy()
    for genes, combos in combo_dicts.items():
        out = pd.Series("Other", index=df.index)
        for combo in combos:
            mask = df["cluster"].str.match(combo_regex[combo])
            assert mask.sum() != 0
            out[mask] = combo
        df[genes] = out
    return df


def make_clusterwise_sc_df(sc_df: pd.DataFrame) -> pd.DataFrame:
    """
    make a mean for each cluster
    """
    out_series: list[pd.Series] = []
    for cluster, sub_df in sc_df.groupby("cluster", observed=True):
        del sub_df["cluster"]
        out_ser = sub_df.mean()
        assert isinstance(out_ser, pd.Series)
        out_ser.name = cluster
        out_series.append(out_ser)
    return pd.DataFrame(out_series)


def plot_genewise_threshold(
    sc_clusters: pd.DataFrame, insitu_expression: dict["str", np.ndarray], scale: float
):
    for gene, expression_array in insitu_expression.items():
        fig, (histax, lineax) = plt.subplots(2, sharex=True)
        assert isinstance(histax, Axes)
        assert isinstance(lineax, Axes)
        histax.set_title(gene)
        histax.hist(expression_array * scale, bins=50)
        lineax.plot(sc_clusters[gene].sort_values(), np.arange(sc_clusters.index.size))


def plot_bars_on_ax(
    clusterwise_sc: pd.DataFrame,
    in_situ_df: pd.DataFrame,
    scale: float,
    cluster_percentages: pd.Series,
    ax: Axes,
):
    in_situ_df = in_situ_df.copy()
    if set(in_situ_df.columns) == set(["Vsx2", "Meis2", "Otx2", "Cell size"]):
        exp_genes = ["Meis2", "Vsx2"]
    else:
        exp_genes = [
            c
            for c in in_situ_df.columns
            if c not in list(BPC_MARKER_THRESHS.keys()) + ["Cell size"]
        ]
    gene1, gene2 = exp_genes
    for gene in exp_genes:
        in_situ_df[f"{gene} on"] = in_situ_df[gene] * scale > gene_threshs[gene]
        clusterwise_sc[f"{gene} on"] = clusterwise_sc[gene] > gene_threshs[gene]
    for df in (in_situ_df, clusterwise_sc):
        for gene, other_gene in zip(exp_genes, exp_genes[::-1]):
            df[f"{gene} alone"] = df[f"{gene} on"] & (~df[f"{other_gene} on"])
        df[f"{gene1} & {gene2}"] = df[f"{gene1} on"] & df[f"{gene2} on"]
        df["Neither"] = (~df[f"{gene1} on"]) & (~df[f"{gene2} on"])
    # make the bar plots
    width = 0.45
    labels = [f"{gene1} & {gene2}", f"{gene1} alone", f"{gene2} alone", "Neither"]
    base_x = np.arange(len(labels))
    # add the in situ data
    offset = width
    bar_height = in_situ_df[labels].mean() * 100
    ax.bar(base_x + offset / 2, bar_height, width, label="In situ", color="k")
    # add the sc data
    bottoms = np.zeros(len(labels))
    for i, (cluster, series) in enumerate(clusterwise_sc.iterrows()):
        brick_height = cluster_percentages[cluster] * 100
        hatch_int, color_int = divmod(i, 7)
        hatch = "/" if hatch_int == 1 else None
        bar_height = series[labels] * brick_height
        ax.bar(
            base_x - offset / 2,
            bar_height,
            width,
            label=cluster,
            bottom=bottoms,
            color=COLORS[color_int],
            hatch=hatch,
        )
        bottoms += bar_height
    # fix labels
    ax.set_xticks(base_x, labels, rotation=70, ha="right", rotation_mode="anchor")


def make_plots():
    in_situ_dfs: list[pd.DataFrame] = []
    for df in get_in_situ_csvs():
        try:
            in_situ_dfs.append(preproc(df))
        except ValueError:
            print(f"{df.columns} has no bpc_marker")
    # get scaling factor between single cell and in_situ
    all_genes = set(np.unique(np.concatenate([d.columns for d in in_situ_dfs])))
    all_genes.discard("Cell size")
    gene_expression_dict: dict[str, np.ndarray] = {}
    for gene in all_genes:
        columns = [l for l in [d.get(gene) for d in in_situ_dfs] if l is not None]
        gene_expression_dict[gene] = np.concatenate(columns)
    # do weighted average of mean fractions
    adata = sc.read_h5ad(HERE / "../data/bpc.h5ad")
    weighted_scales = [
        in_situ_exp.size * adata[:, gene].X.mean() / in_situ_exp.mean()
        for gene, in_situ_exp in gene_expression_dict.items()
    ]
    total_counts = sum(a.size for a in gene_expression_dict.values())
    print(f"{total_counts=}")
    scale = sum(weighted_scales) / total_counts

    # plot the dfs
    sc_df = adata[:, list(all_genes)].to_df()
    sc_df["cluster"] = adata.obs["cluster"]
    clusterwise_sc = make_clusterwise_sc_df(sc_df)
    plot_genewise_threshold(clusterwise_sc, gene_expression_dict, scale)
    cluster_percentages = sc_df.cluster.value_counts() / sc_df.index.size
    with plt.style.context(HERE / "paper.mplstyle"):
        fig, axs = plt.subplots(1, 5, sharey=True, layout="constrained")
    fig.set_size_inches((16, 6.3))
    axs[0].set_ylabel("Percent of bipolar cells")
    for ax, in_situ_df in zip(axs, in_situ_dfs):
        plot_bars_on_ax(clusterwise_sc, in_situ_df, scale, cluster_percentages, ax)
    fig.legend(
        *axs[0].get_legend_handles_labels(),
        bbox_to_anchor=(0.5, -0.1),
        loc="lower center",
        borderaxespad=0,
        fontsize=15,
        ncol=8,
    )
    sns.despine(fig)
    fig.savefig(HERE / "../imgs/bpcs_quant.svg")
    # avg_sc_df = combine_groups(average_sc_rnaseq(sc_df, 100, 10))
    # for df in in_situ_dfs:
    # fig, ax = plt.subplots()
    # plot_in_situ_on_axs(df, avg_sc_df, ax, scale, debug_plot=False)


def plot_otxvsx():
    (raw_df,) = get_in_situ_csvs(
        {
            HERE
            / "../quantified_image_data/Well1_405wga_488Vsx2_546Meis2_647Otx2-1.csv": (
                "Vsx2",
                "Meis2",
                "Otx2",
            ),
        }
    )
    in_situ_df = preproc(raw_df, False)
    exp_genes = "Vsx2", "Otx2"
    gene1, gene2 = exp_genes
    for gene in exp_genes:
        in_situ_df[f"{gene} on"] = in_situ_df[gene] > bpc_marker_threshs[gene]
    for gene, other_gene in zip(exp_genes, exp_genes[::-1]):
        in_situ_df[f"{gene} alone"] = in_situ_df[f"{gene} on"] & (
            ~in_situ_df[f"{other_gene} on"]
        )
    in_situ_df[f"{gene1} & {gene2}"] = (
        in_situ_df[f"{gene1} on"] & in_situ_df[f"{gene2} on"]
    )
    in_situ_df["Neither"] = (~in_situ_df[f"{gene1} on"]) & (~in_situ_df[f"{gene2} on"])
    in_situ_df.to_csv(HERE / "../quantified_image_data/bpc_marker.csv")
    # make the bar plots
    width = 0.9
    labels = [f"{gene1} & {gene2}", f"{gene1} alone", f"{gene2} alone", "Neither"]
    base_x = np.arange(len(labels))
    # add the in situ data
    bar_height = in_situ_df[labels].mean() * 100
    with plt.style.context("paper.mplstyle"):
        fig, ax = plt.subplots(layout="constrained")
    fig.set_size_inches((3.4, 6.3))
    ax.bar(base_x, bar_height, width, label="In situ", color="k")
    ax.set_xticks(base_x, labels, rotation=70, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Percent of cells in INL")
    sns.despine(fig)
    fig.savefig(HERE / "../data/all_cells.svg")


if __name__ == "__main__":
    plot_otxvsx()
    make_plots()
