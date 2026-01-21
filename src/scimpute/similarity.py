import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import os
from pathlib import Path
from scimpute.utils import intersection, checkformat
from scimpute.io import read_matrix, read_cell_identities

def similarity_matrix(
    x,
    y,
    clusters,
    cell_name_column_idx,
    outdir,
    metric = "cosine_similarity",
    k_neighbors = 25,
    consider_clusters = True
):
    """
    Calculates similarity matrix between the cells of two single cell transcriptomics datasets.

    Parameters
    ----------
    x : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to file containing the gene expression matrix to perform gene expression imputation for (e.g. from a spatial experiment); genes as rows, cells/spots as columns
    y : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to file containing the comprehensive gene expression matrix to perform gene expression imputation from (e.g. from a single cell/nucleus RNA sequencing experiment); genes as rows, cells/spots as columns
    clusters : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to file containing the identities for all cells from both datasets
    cell_name_column_idx : int, required
        index of the column containing the cell names (not the cell type clusters)
    outdir : str or pathlib.Path, required
        Path to a directory where the similarity matrix determination output should be written on the disk.
    metric : {"cosine_similarity"}
        similarity metric (currently only supporting cosine similarity); default is 'cosine_similarity'
    k_neighbors : int, optional
        number of k nearest neighbors to find between datasets; default is '25'
    consider_clusters : bool, optional
        flag indicating whether to only identify similar cells between datasets within the same annotated cluster; default is 'True'
    """

    xformat = checkformat(x)
    yformat = checkformat(y)
    clustersformat = checkformat(clusters)

    if xformat == "path":
        x = read_matrix(x)
    if yformat == "path":
        y = read_matrix(y)
    if clustersformat == "path":
        clusters = read_cell_identities(clusters, cell_name_column_idx)

    Path(outdir).mkdir(parents=True, exist_ok=True)

    y_red = y[y.columns.intersection(x.columns)]
    x = x[y_red.columns]
    y_red = y_red[x.columns]

    unique_clusters = clusters["id"].unique()
    result_tuples = []


    for cluster_id in unique_clusters:
        if consider_clusters == True:
            x_cluster = x[clusters["id"] == cluster_id]
            y_cluster = y_red[clusters["id"] == cluster_id]
        else:
            x_cluster = x
            y_cluster = y_red
        
        if metric == "cosine_similarity":
            cosine_sim_matrix = cosine_similarity(x_cluster, y_cluster)

        distance_matrix_cluster_df = pd.DataFrame(cosine_sim_matrix, index=x_cluster.index, columns=y_cluster.index)

        for row_name, row in distance_matrix_cluster_df.iterrows():
            top_indices = (-row).argsort()[:k_neighbors] # from, to, -1 indicates descending order, getting the largest entries here
            top_values = row.iloc[top_indices]
            top_columns = top_values.index
            result_tuples.extend([(cluster_id, row_name, col_name, value+1) for col_name, value in zip(top_columns, top_values)])
    result_df = pd.DataFrame(result_tuples, columns=["cluster", "x", "y", "distance"])
    result_df.to_csv(os.path.join(outdir, "sim_matrix.tsv"), sep="\t", index=False)
    return result_df