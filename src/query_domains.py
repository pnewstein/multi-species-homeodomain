from pathlib import Path
import json
import time
from io import StringIO
from typing import Iterable
from urllib.parse import quote

import anndata as ad
import numpy as np
import pandas as pd
import requests

try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path()

RATE_LIMIT = 2

homeodomain_smart_code = "SM00389"
mouse_dset = "mmusculus_gene_ensembl"
fly_dset = "dmelanogaster_gene_ensembl"
biomart_template = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >

<Dataset name = "{dset}" interface = "default" >
<Filter name = "ensembl_gene_id" value = "{id}"/>
<Attribute name = "ensembl_gene_id" />
<Attribute name = "hmmpanther" />
<Attribute name = "smart" />
</Dataset>
</Query>
"""


def get_hdtf_status(ids: Iterable[str], fly=False) -> dict[str, tuple[bool, bool]]:
    """
    Uses biomart to query the smart annotations for this gene returnes
    a dict mapping each id to whether it has a homeodomain
    """
    cache_path = HERE / "../data/id_to_hdtf_cache.json"
    try:
        out: dict[str, tuple[bool, bool]] = json.loads(cache_path.read_text())
    except FileNotFoundError:
        out = {}
    ids_to_check = set(ids) - set(out.keys())
    if len(ids_to_check) != 0:
        panther_hdtf_families = np.array(
            pd.read_csv(HERE / "../data/panther_hdtf_families.txt", header=None)[0]
        )
        url = "https://www.ensembl.org/biomart/martservice"
        headers = {"Content-Type": "application/xml"}
        dset = fly_dset if fly else mouse_dset
        for chunk in np.array_split(list(ids_to_check), len(ids_to_check) // 500 + 1):
            xml = biomart_template.format(id=",".join(chunk), dset=dset)
            time.sleep(RATE_LIMIT)
            response = requests.post(url, data={"query": xml, "headers": headers})
            while response.status_code != 200:
                time.sleep(60)
                response = requests.post(url, data={"query": xml, "headers": headers})
            string_io = StringIO(response.text)
            string_io.seek(0)
            response_df = pd.read_table(string_io, delimiter="\t")
            response_df["panther_hdtf"] = response_df["PANTHER ID"].isin(
                panther_hdtf_families
            )
            response_df["smart_hdtf"] = (
                response_df["SMART ID"] == homeodomain_smart_code
            )
            for gid in chunk:
                gene_annotations = response_df.loc[response_df["Gene stable ID"] == gid]
                if len(gene_annotations) == 0:
                    out[gid] = (False, False)
                else:
                    out[gid] = (
                        bool(gene_annotations["smart_hdtf"].max()),
                        bool(gene_annotations["panther_hdtf"].max()),
                    )
        cache_path.write_text(json.dumps(out))
    # return only those ids your queried
    return {id: out[id] for id in ids}


def get_expressed_genes(
    adata: ad.AnnData, transcript_thresh: float, obs_name: str
) -> tuple[set[str], dict[str, np.ndarray]]:
    """
    gets all genes expressed above a certain thresh in a category in the obs table
    """
    all_expressed_genes: set[str] = set()
    gene_expression: dict[str, np.ndarray] = {}
    for cluster, sub_df in adata.obs.groupby(obs_name):
        assert isinstance(cluster, str)
        sub_adata = adata[sub_df.index, :]
        mean_expression_array_like = sub_adata[sub_df.index, :].X.mean(axis=0)
        if isinstance(mean_expression_array_like, np.matrix):
            mean_expression_array = np.array(mean_expression_array_like).flatten()
        else:
            mean_expression_array = mean_expression_array_like
        is_expressed = mean_expression_array > transcript_thresh
        expressed_genes = adata.var_names[is_expressed]
        assert isinstance(expressed_genes, pd.Index)
        all_expressed_genes.update(expressed_genes)
        gene_expression[cluster] = np.array(expressed_genes.values)
    return all_expressed_genes, gene_expression


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


def search_extern_name(symbols: list[str]) -> tuple[dict[str, str], list[str]]:
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


def get_ensambl_ids(expressed_genes: pd.Series) -> tuple[dict[str, str], list[str]]:
    """
    gets the ensambl ides of each gene
    """
    known_unknowns_path = Path(HERE / "../data/known_unknowns.json")
    try:
        unknowns: list[str] = json.loads(known_unknowns_path.read_text())
    except FileNotFoundError:
        unknowns = []
    cache_path = Path(HERE / "../data/gene_ids.json")
    try:
        name_to_id_series = pd.Series(json.loads(cache_path.read_text()))
    except FileNotFoundError:
        name_to_id_series = pd.Series()
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
            name_to_id_series = pd.concat(
                [name_to_id_series, pd.Series(id_dict | searched_id_dict)]
            )
            cache_path.write_text(json.dumps(name_to_id_series.to_dict()))
            known_unknowns_path.write_text(json.dumps(unknowns))
    return name_to_id_series, unknowns


def get_expressed_hdtf(
    adata: ad.AnnData, transcript_thresh: float, obs_name: str, fly=False
) -> dict[str, np.ndarray[str]]:
    genes, cluster_expression = get_expressed_genes(adata, transcript_thresh, obs_name)
    if fly:
        id_series = pd.read_csv(HERE / "../data/gene_ids.csv", index_col=0)["gene_ids"]
        unknowns = []
    else:
        # use ensambl to seach the ensambl id from the name
        id_series, unknowns = get_ensambl_ids(pd.Series(list(genes)))
        print(unknowns)
    # exclude unknowns
    genes = genes - set(unknowns)
    ids = id_series[list(genes)]
    assert isinstance(ids, pd.Series)
    hdtf_status_dict = get_hdtf_status(ids, fly)
    hdtf_df = pd.DataFrame(hdtf_status_dict).sum().astype(bool)
    # filter those genes which are not HDTFs
    cluster_hdtf_expression = {}
    for cluster, cluster_expressed in cluster_expression.items():
        known_cluster_expression = list(set(cluster_expressed) - set(unknowns))
        cluster_ids = id_series[known_cluster_expression]
        name_to_hdtf = cluster_ids.map(hdtf_df)
        cluster_hdtf_expression[cluster] = np.array(
            name_to_hdtf.loc[name_to_hdtf].index
        )
    return cluster_hdtf_expression
