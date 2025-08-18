"""
Makes all of the filtered h5ad for dotplots
"""

from pathlib import Path
import pandas as pd
import json
import time
from io import StringIO
from itertools import chain
from urllib.parse import quote
import re
from typing import Sequence

import scanpy as sc
import requests
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import umap

try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()

RATE_LIMIT = 2
major_classes = ("BC", "Cone", "Rod", "HC", "AC", "")
homeodomain_smart_code = "SM00389"
biomart_template = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >

<Dataset name = "mmusculus_gene_ensembl" interface = "default" >
<Filter name = "ensembl_gene_id" value = "{id}"/>
<Attribute name = "ensembl_gene_id" />
<Attribute name = "ensembl_gene_id_version" />
<Attribute name = "smart" />
</Dataset>
</Query>
"""


def get_ensamble_gene_info(symbols: list[str]) -> dict[str, str]:
    """
    maps all symbols to their ids
    """
    print(f"getting gene info for {len(symbols)} genes")
    species = "Mus_musculus"
    url = f"https://rest.ensembl.org/lookup/symbol/{species}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    time.sleep(RATE_LIMIT)
    responce = requests.post(
        url, headers=headers, data=json.dumps({"symbols": symbols})
    )
    return {name: table["id"] for name, table in responce.json().items()}

def search_extern_name(symbols: list[str]) -> dict[str, str]:
    """
    gets best guess for external name for gene must do one request per gene
    """
    species = "Mus_musculus"
    headers = {"Content-Type": "application/json"}
    out: dict[str, str] = {}
    unknown = []
    for gene_name in symbols:
        time.sleep(RATE_LIMIT)
        url = f"https://rest.ensembl.org/xrefs/symbol/{species}/{quote(gene_name)}?"
        response = requests.get(url, headers=headers)
        genes = [e for e in response.json() if e["type"] == "gene"]
        if len(genes) != 1:
            print(len(genes), gene_name)
            unknown.append(gene_name)
        else:
            out[gene_name] = genes[0]["id"]
    return out, unknown
    

def get_hdtf_status(ids: list[str]) -> dict[str, bool]:
    """
    Uses biomart to query the smart annotations for this gene returnes
    a dict mapping each id to whether it has a homeodomain
    """
    print(f"getting hdtf status for {len(ids)} genes")
    time.sleep(RATE_LIMIT)
    url = "https://www.ensembl.org/biomart/martservice"
    headers = {"Content-Type": "application/xml"}
    xml = biomart_template.format(id=",".join(ids))
    response = requests.post(url, data={"query": xml, "headers": headers})
    while response.status_code != 200:
        time.sleep(60)
        response = requests.post(url, data={"query": xml, "headers": headers})
    string_io = StringIO(response.text)
    string_io.seek(0)
    response_df = pd.read_table(string_io, delimiter="\t")
    out = {s: False for s in ids}
    if "SMART ID" not in response_df.columns:
        print(response_df)
        return out
    for _, row in response_df.iterrows():
        if row["SMART ID"] == homeodomain_smart_code:
            out[row.at["Gene stable ID"]] = True
    return out

def get_expressed_genes(adata: sc.AnnData, transcript_thresh: float, obs_name: str) -> set[str]:
    """
    gets all genes expressed above a certain thresh in a category in the obs table
    """
    all_expressed_genes: set[str] = set()
    for _, sub_df in adata.obs.groupby(obs_name):
        sub_adata = adata[sub_df.index, :]
        mean_expression_array_like = sub_adata[sub_df.index, :].X.mean(axis=0)
        if isinstance(mean_expression_array_like, np.matrix):
            mean_expression_array = np.array(mean_expression_array_like).flatten()
        else:
            mean_expression_array = mean_expression_array_like
        is_expressed = mean_expression_array > transcript_thresh
        expressed_genes = adata.var_names[is_expressed]
        all_expressed_genes.update(expressed_genes)
    return all_expressed_genes

genes = ['Skp1a', 'E130218I03Rik', 'Atp5e', 'Atp5l', 'Hist3h2ba', 'Atp5mpl', 'Atp5f1', 'C130071C03Rik', 'Atp5j', 'Mpp6', 'Atp5j2', 'Atp5g2', 'Gm45716', 'Hist3h2a', 'Fam172a', 'Fam189a1', 'Sept4', 'Hist1h2b c', 'Gm48693', 'A730046J19Rik', 'Atp5c1', 'Kirrel', 'A530058N18Rik', 'Atp5o', 'Atp5md', 'Qk', 'Atp5g3', 'Nars', 'Ndufb1-ps', 'Gm20754', 'Hprt', 'Gm4258', 'Atp5g1', 'H2afz', 'H2afy', 'Atp5h', 'Pnmal2', 'March1 ', 'AC149090.1', 'Gm45323', 'March3', 'Atp5b', 'Atp5d', 'Atp5a1', 'Atp5k', 'Fam155a', 'Sept7', 'Atpif1', 'Gm42418']

def get_ensambl_ids(expressed_genes: pd.Series) -> dict[str, str]:
    """
    gets the ensambl ides of each gene
    """
    known_unknowns_path = Path(HERE / "../data/known_unknowns.json")
    # known_unknowns_path.write_text(json.dumps(unknowns))
    try:
        unknowns: list[str] = json.loads(known_unknowns_path.read_text())
    except FileNotFoundError:
        unknowns = []
    cache_path = Path(HERE / "../data/gene_ids.json")
    try:
        name_to_id_series = pd.Series(json.loads(cache_path.read_text()))
    except FileNotFoundError:
        name_to_id_series = pd.Series()
    # expressed_genes = pd.Series(gene_list)
    all_unknown_genes = expressed_genes[~expressed_genes.isin(name_to_id_series.index)]
    # filter those previously known to be unknown
    unknown_genes = all_unknown_genes.loc[~all_unknown_genes.isin(unknowns)]
    # enure name_to_hdtf_series contains unknown_genes
    if len(unknown_genes) != 0:
        for chunk in np.array_split(unknown_genes, len(unknown_genes) // 500 + 1):
            id_dict = get_ensamble_gene_info(chunk.to_list())
            # individualy search each unknown gene to try to determine id
            missing_genes = chunk.loc[~chunk.isin(id_dict.keys())]
            searched_id_dict, unknown = search_extern_name(missing_genes)
            unknowns.extend(unknown)
            name_to_id_series = pd.concat([name_to_id_series, pd.Series(id_dict | searched_id_dict)])
            cache_path.write_text(json.dumps(name_to_id_series.to_dict()))
    return name_to_id_series, unknowns

def get_hdtf_lists(gene_sets: dict[str, set[str]]) -> dict[list[str]]:
    """
    gets all of the HDTFS from all the genes
    """
    # denest the list to reduce api calls
    expressed_genes = pd.Series(list(set(chain.from_iterable(gene_sets.values()))))
    name_to_id_series, unknowns = get_ensambl_ids(expressed_genes)
    hdtf_df = compare_hdtf_lists()
    hdtf_df = hdtf_df.drop(["superfamily", "interpro"], axis=1)
    hdtf_df = hdtf_df.loc[hdtf_df.sum(axis=1) > 0]
    name_to_hdtf_series = pd.Series(np.array(expressed_genes.isin(hdtf_df.index)), index=expressed_genes)
    print(expressed_genes[expressed_genes.isin(unknowns)])
    hdtf_lists: dict[str, list[str]] = {}
    for cell_class, gene_set in gene_sets.items():
        hdtf_series = name_to_hdtf_series[list(gene_set.difference(unknowns))]
        hdtf_list = hdtf_series[hdtf_series].index.to_list()
        hdtf_lists[cell_class] = hdtf_list
    return hdtf_lists

def get_smart_hdtfs(gene_sets: dict[str, set[str]]) -> dict[list[str]]:
    """
    gets all of the HDTFS from all the genes
    """
    # denest the list to reduce api calls
    expressed_genes = pd.Series(list(set(chain.from_iterable(gene_sets.values()))))
    # try reading from cache
    known_unknowns_path = Path(HERE / "../data/known_unknowns.json")
    # known_unknowns_path.write_text(json.dumps(unknowns))
    try:
        unknowns: list[str] = json.loads(known_unknowns_path.read_text())
    except FileNotFoundError:
        unknowns = []

    cache_path = Path(HERE / "../data/hdtf_genes.json")
    try:
        name_to_hdtf_series = pd.Series(json.loads(cache_path.read_text()))
    except FileNotFoundError:
        name_to_hdtf_series = pd.Series().astype(bool)
    all_unknown_genes = expressed_genes[~expressed_genes.isin(name_to_hdtf_series.index)]
    # filter those previously known to be unknown
    unknown_genes = all_unknown_genes.loc[~all_unknown_genes.isin(unknowns)]
    # enure name_to_hdtf_series contains unknown_genes
    if len(unknown_genes) != 0:
        for chunk in np.array_split(unknown_genes, len(unknown_genes) // 200 + 1):
            id_series = pd.Series(get_ensamble_gene_info(chunk.to_list()))
            # individualy search each unknown gene to try to determine id
            missing_genes = chunk.loc[~chunk.isin(id_series.index)]
            searched_id_dict, unknown = search_extern_name(missing_genes)
            unknowns.extend(unknown)
            id_series = pd.concat([id_series, pd.Series(searched_id_dict)])
            hdtf_series = pd.Series(get_hdtf_status(list(id_series)))
            new_name_to_hdtf_series = hdtf_series.reindex(id_series.values).set_axis(id_series.index)
            name_to_hdtf_series = pd.concat([new_name_to_hdtf_series, name_to_hdtf_series])
            # cache the result
            cache_path.write_text(json.dumps(name_to_hdtf_series.to_dict()))
    hdtf_lists: dict[str, list[str]] = {}
    for cell_class, gene_set in gene_sets.items():
        hdtf_series = name_to_hdtf_series[list(gene_set.difference(unknowns))]
        hdtf_list = hdtf_series[hdtf_series].index.to_list()
        hdtf_lists[cell_class] = hdtf_list
    return hdtf_lists


def plot_dotplot(adata: sc.AnnData, obs_name: str, thresh: float, cmap: str, order: Sequence[str] | None =None, use_umap=False, fly=False):
    cluster_hdtf_expression = get_expressed_hdtf(adata, thresh, obs_name, fly)
    genes = list(set(chain.from_iterable(cluster_hdtf_expression.values())))
    cluster_names = np.unique(adata.obs[obs_name])
    expression_df = adata[:, np.array(genes)].to_df()
    expression_df[obs_name] = adata.obs[obs_name]
    clusterwise_expression = pd.DataFrame(np.nan, index=cluster_names, columns=genes)
    for cluster, sub_df in expression_df.groupby(obs_name, observed=True):
        if cluster == "other": continue
        clusterwise_expression.loc[cluster, :] = sub_df.loc[:, genes].mean()
    if order is None:
        number_above_thresh = (clusterwise_expression > thresh).sum()
        sum_expression = clusterwise_expression.sum(axis=0)
        hdtfs_in_order = pd.DataFrame([number_above_thresh, sum_expression]).T.sort_values([0, 1]).index
    else:
        hdtfs_in_order = list(order)
    clusters_in_order = np.sort(clusterwise_expression.index).tolist()
    if use_umap:
        reduced = umap.UMAP(n_components=1, random_state=100).fit_transform(clusterwise_expression).flatten()
        clusters_in_order = (
            pd.Series(reduced, index=cluster_names)
            .sort_values()
            .index
        )
    ax: Axes
    fig: Figure
    fig, ax = plt.subplots(constrained_layout=True) # type: ignore
    edge_color = ax.spines["left"].get_edgecolor()
    dp = sc.pl.DotPlot(adata, var_names=hdtfs_in_order, groupby=obs_name, categories_order=clusters_in_order, ax=ax)
    dp.style(cmap=cmap, dot_edge_color=edge_color)
    dp.swap_axes()
    axs = dp.show(return_axes=True)
    sns.despine(fig)
    return fig, axs

def get_gene_counts(gene: str, majorclass: str, adata: sc.AnnData):
    sub_adata = adata[adata.obs["majorclass"] == majorclass, gene]
    out = pd.Series({channel: sub_adata[sub_df.index].X.mean() for channel, sub_df in sub_adata.obs.groupby("celltype", observed=True)})
    return out

def get_hdtf_from_fasta(name_to_id_series: pd.Series) -> pd.Series:
    """
    takes the name to id series mapping names to ids
    reads the fasta file for gene ids
    returns a boolian listing of which genes are hdtfs
    """
    regex = re.compile(r"^>.*gene:(ENSMUSG\d+) ")
    text = Path(HERE / "../data/Mus_musculus_46689.fasta").read_text()
    gene_ids: set[str] = set()
    matches = re.finditer(r"^>.*gene:(ENSMUSG\d+) ", text, re.MULTILINE)
    for re_match in matches:
        (gene_id,) = re_match.groups()
        gene_ids.add(gene_id)
    out = name_to_id_series.isin(gene_ids)
    return out


def get_ensambl_ids_from_json_list(json_dict: dict[str, list[dict[str, str]]]):
    return [e["to"].split(".")[0] for e in json_dict["results"]]

def compare_hdtf_lists():
    cache_path = Path(HERE / "../data/gene_ids.json")
    name_to_id_series = pd.Series(json.loads(cache_path.read_text()))
    interpro_list = json.loads(Path("interpro.json").read_text())
    interpro_ids = get_ensambl_ids_from_json_list(interpro_list)
    interpro_names = name_to_id_series.loc[name_to_id_series.isin(interpro_ids)]
    panther_list = json.loads(Path("panther_ids.json").read_text())
    panther_ids = get_ensambl_ids_from_json_list(panther_list)
    panther_names = name_to_id_series.loc[name_to_id_series.isin(panther_ids)]
    superfamily_names = name_to_id_series.loc[get_hdtf_from_fasta(name_to_id_series)]
    name_to_hdtf_series = pd.Series(json.loads(Path(HERE / "../data/hdtf_genes.json").read_text()))
    smart_df = pd.DataFrame({"id": name_to_id_series, "smart": name_to_hdtf_series})
    smart_names = smart_df.loc[smart_df["smart"], "id"]
    full_set = pd.Index(np.unique(np.concatenate((smart_names.index, superfamily_names.index, panther_names.index, interpro_names.index))))
    out_df = pd.DataFrame(index=full_set, data={"interpro": full_set.isin(interpro_names.index), "panther": full_set.isin(panther_names.index), "superfamily": full_set.isin(superfamily_names.index), "smart": full_set.isin(smart_names.index)})
    # out_df = out_df.drop(["interpro"], axis=1)
    #out_df = out_df.loc[[np.std(r[1]) != 0 for r in out_df.iterrows()]]
    # out_df = out_df.loc[[np.std(r[1]) != 0 for r in out_df.iterrows()]]
    # out_df = out_df.loc[out_df.index.isin(expressed_genes)]

    return out_df



cells_we_looked_at = ["Irx5", "Lhx3", "Sebox", "Zeb2", "Six3", "Meis2", "Vsx1", "Zfhx4", "Isl1", "Lhx4", "Vsx2", "Otx2"]
def main():
    adata = sc.read_h5ad(HERE / "../data/MRCA_full_normcounts.h5ad")
    adata = adata[~adata.obs["celltype"].isin(["RGC", "AC", "BC"])]
    adata.obs["majorclass"] =  adata.obs["majorclass"].astype(str)
    adata.obs.loc[adata.obs["majorclass"].isin(["Rod", "Cone"]), "majorclass"] = "PR"
    adata.obs["majorclass"] = adata.obs["majorclass"].astype("category")
    series = adata.obs.groupby(["majorclass", "celltype"])["age"].count()
    counts = series.loc[series != 0]
    print(counts)
    thresh = .35
    classes = ["AC", "RGC", "BC", "PR", "HC"]
    gene_sets: dict[str, set[str]] = {}
    for cell_class in classes:
        adata_in_class = adata[adata.obs["majorclass"] == cell_class]
        gene_sets[cell_class] = get_expressed_genes(adata_in_class, thresh, "celltype")
    hdtf_lists = get_hdtf_lists(gene_sets)
    print(np.isin(cells_we_looked_at, hdtf_lists["BC"]))
    print(len(hdtf_lists["BC"]))
    size_dict = {
            "AC": (16.61, 7.49),
            "RGC": (12.5, 7.77),
            "BC": (12.85, 7.65),
            "PR": (2.02, 3.83),
            "HC": (1.86, 3.2)
    }
    figures = []
    for cell_class, hdtf_list in hdtf_lists.items():
        size = size_dict[cell_class]
        adata_in_class = adata[adata.obs["majorclass"] == cell_class]
        use_pca = False if cell_class == "PR" else True
        with plt.style.context("style.mplstyle"):
            print(plt.rcParams['font.family'])
            fig, axs = plot_dotplot(adata_in_class, "celltype", hdtf_list, thresh=0.2, cmap="plasma_r", use_pca=use_pca)
        fig.set_size_inches(size)
        fig.savefig(HERE / f"imgs/{cell_class}.svg")
        figures.append(fig)



