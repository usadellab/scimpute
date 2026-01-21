import pandas as pd
import os
from pathlib import Path
from math import ceil
from scimpute.utils import intersection, checkformat
from scimpute.io import read_matrix, read_matrix_basic


def impute_expression(
    x,
    y,
    similarity_matrix,
    outdir,
    chunk_size = 1000,
):
    """
    Runs gene expression imputation chunk-wise.

    Parameters
    ----------
    x : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to file containing the gene expression matrix to perform gene expression imputation for (e.g. from a spatial experiment); genes as rows, cells/spots as columns
    y : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to file containing the comprehensive gene expression matrix to perform gene expression imputation from (e.g. from a single cell/nucleus RNA sequencing experiment); genes as rows, cells/spots as columns
    similarity_matrix : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to file containing cluster identity table, cell names from the dataset to impute gene expression for, cell names from the dataset to impute gene expression from, cell similarity/distance value
    outdir : str or pathlib.Path, required
        location of output directory
    chunk_size : int, optional
        number of cells per chunk
    """

    xformat = checkformat(x)
    yformat = checkformat(y)
    similarity_matrixformat = checkformat(similarity_matrix)

    if xformat == "path":
        x = read_matrix(x)
    if yformat == "path":
        y = read_matrix(y)
    if similarity_matrixformat == "path":
        similarity_matrix = read_matrix_basic(similarity_matrix)

    Path(os.path.join(outdir, "imputations")).mkdir(parents=True, exist_ok=True)

    imputation_cell_idx_start = 0
    imputation_cell_idx_stop = imputation_cell_idx_start + chunk_size
    nIterations = ceil(len(x.index) / chunk_size)


    for i in range(nIterations):
        imputed_results = []
        if imputation_cell_idx_stop > len(x.index):
            imputation_cell_idx_stop = len(x.index)
        print(f"imputation in range of: {imputation_cell_idx_start} to {imputation_cell_idx_stop}")

        for target_cell in x.index[imputation_cell_idx_start:imputation_cell_idx_stop]:
            distances_for_target = similarity_matrix[similarity_matrix["x"] == target_cell]
            neighboring_cells = distances_for_target["y"]
            distances = distances_for_target["distance"]
        
            weights = distances  # Inverse distance as weights (you may choose a different weighting scheme)
            weights = weights.reset_index(drop=True)
        
        
            # Get neighboring gene expressions
            neighboring_gene_expression = y.loc[neighboring_cells]
            neighboring_gene_expression = neighboring_gene_expression.reset_index(drop=True)
        
            # Calculate imputed expression for each gene
            imputed_expression = round((neighboring_gene_expression.T * weights).sum(axis=1) / weights.sum(), 2)
        
            imputed_results.append({"cell": target_cell, **imputed_expression.to_dict()})

        columns = ["cell"] + list(y.columns)
        imputed_results_df = pd.DataFrame(imputed_results, columns=columns)
        imputed_results_df.to_csv(os.path.join(outdir, "imputations", f"CPM_imputation_{i:03}_{imputation_cell_idx_start}_{imputation_cell_idx_stop}.tsv"), sep="\t")

        imputation_cell_idx_start = imputation_cell_idx_stop
        imputation_cell_idx_stop = chunk_size * (i+2)
    
    return imputed_results_df