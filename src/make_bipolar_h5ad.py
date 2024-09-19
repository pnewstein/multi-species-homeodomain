"""
this code creates an H5ad from expression matrix and metadata
"""

from pathlib import Path
import re

import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix

try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()


def main():
    metadata = pd.read_csv(
        HERE / "../seq-data/clust_retinal_bipolar.txt", index_col=0, header=1, sep="\t"
    )
    metadata.columns = ["cluster", "subcluster"]
    assert np.all(metadata["cluster"] == metadata["subcluster"])
    matrix = pd.read_csv(
        HERE / "../seq-data/exp_matrix.txt", index_col=0, header=0, sep="\t"
    ).T
    adata = sc.AnnData(matrix)
    adata.X = csr_matrix(adata.X)
    # get rid of other celltypes
    adata.obs["cluster"] = metadata["cluster"].astype("category")
    is_bpc_mask = np.logical_not(
        adata.obs["cluster"].isin(
            [
                "Doublets/Contaminants",
                "AC (Amacrine cell)",
                "Cone Photoreceptors",
                "MG (Mueller Glia)",
                "Rod Photoreceptors",
            ]
        )
    )
    filtered_adata = adata[is_bpc_mask]
    # strip out parantheticals
    filtered_adata.obs["cluster"] = filtered_adata.obs["cluster"].cat.rename_categories(
        [
            re.sub(r"\(.+\)", "", cat).strip()
            for cat in filtered_adata.obs["cluster"].cat.categories
        ]
    )
    # do simple umap
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    filtered_adata.write_h5ad(HERE / "../data/bpc.h5ad")


if __name__ == "__main__":
    main()
