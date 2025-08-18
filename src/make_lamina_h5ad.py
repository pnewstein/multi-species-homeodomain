"""
Creates a normalized H5AD file from 10x cellranger output
"""

from pathlib import Path
import itertools

import scanpy as sc
import anndata as ad
import pandas as pd
from samalg import SAM
import numpy as np

def main():
    data_path = Path("/mnt/f/peter")
    adata_r1 = sc.read_10x_mtx(data_path/"adult_lamina_scrnaseq/Count/R1/outs/raw_feature_bc_matrix", cache=False)
    adata_r2 = sc.read_10x_mtx(data_path/"adult_lamina_scrnaseq/Count/R2/outs/raw_feature_bc_matrix", cache=False)
    adata_r1.var["gene_ids"].to_csv("gene_ids.csv")
    adata = ad.concat({"r1": adata_r1, "r2": adata_r2}, join="outer", label="dataset_name")
    adata.obs_names_make_unique()
    adata.var['mt'] = adata.var_names.str.startswith('mt:')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    adata.layers["counts"] = adata.X.copy()
    sam = SAM(counts=adata)
    sam.preprocess_data()
    sam.run()
    DIFF_LAMINA_GENES = ["ap","pdm3","bsh", "zfh1", "erm", "chp", "bab2", "svp",
                         "scro", "repo", "nSyb", "VGlut"]
    sc.pl.umap(sam.adata, color=DIFF_LAMINA_GENES)
    sam.clustering(param=1)
    sc.pl.umap(sam.adata, color=["leiden_clusters"], legend_loc="on data", add_outline=True)
    sc.pl.dotplot(sam.adata, var_names=DIFF_LAMINA_GENES, groupby="leiden_clusters")
    cluster_identiy_map = {
        "L1": (19, ), #svp+ zfh1+ VGlut+
        "L2": (1,), #bab2
        "L3": (4, ), #zfh1+ erm+
        "L4": (11,), #Bsh+ Ap+
        "L5": (6, ), #Bsh+ Pdm3+
    }
    clusters_to_keep = list(itertools.chain.from_iterable(cluster_identiy_map.values()))
    adata.obs["leiden_clusters"] = sam.adata.obs["leiden_clusters"]
    l_adata = adata[sam.adata.obs["leiden_clusters"].isin(clusters_to_keep)]
    inverted = pd.Series({v: k for k, vals in cluster_identiy_map.items() for v in vals})
    l_adata.obs["category"] = np.array(inverted[l_adata.obs["leiden_clusters"]])
    sam = SAM(counts=l_adata)
    sam.preprocess_data()
    sam.run()
    sam.clustering(param=.1)
    sc.pl.umap(sam.adata, color="leiden_clusters")
    sc.pl.umap(sam.adata, color=DIFF_LAMINA_GENES)
    sam.adata.write_h5ad("chundi_lamina.h5ad")


def test():
    ca = sc.read_h5ad("../chundi_lamina.h5ad")
    # sam.adata.obs["is_old"] = sam.adata.obs_names.isin(ca.obs_names)
    # sc.pl.umap(sam.adata, color=["is_old"])
