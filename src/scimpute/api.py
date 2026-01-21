from scimpute.similarity import similarity_matrix
from scimpute.imputation import impute_expression
from scimpute.merge import merge_imputation_chunks
from scimpute.io import read_inputs
from scimpute.validation import validate_results
import os
from pathlib import Path


def expression_imputation(
    matrix_to_impute_for,
    matrix_to_impute_from,
    cell_identities,
    cell_name_column_idx,
    outdir = "imputation_output",
    metric = "cosine_similarity",
    k_neighbors = 25,
    consider_clusters = True,
    save_chunks = True,
    chunk_size = 1000
):
    """
    This tool will perform gene expression imputation for a dataset with only limited gene expression information (such as from a spatial dataset) based on a dataset with comprehensive gene expression information (such as from single cell/nucleus RNA sequencing).

    Parameters
    ----------
    matrix_to_impute_for : str or pathlib.Path, required
        tab-separated gene expression matrix to perform gene expression imputation for (e.g. from a spatial experiment); genes as rows, cells/spots as columns
    matrix_to_impute_from : str or pathlib.Path, required
        tab-separated comprehensive gene expression matrix to perform gene expression imputation from (e.g. from a single cell/nucleus RNA sequencing experiment); genes as rows, cells/spots as columns
    cell_identities : str or pathlib.Path, required
        tab-separated table containing identities for all cells from both datasets
    cell_name_column_idx : int, required
        index of the column containing the cell names (not the cell type clusters)
    outdir : str or pathlib.Path, optional
        location of output directory
    metric : {"cosine_similarity"}, optional
        similarity metric (currently only supporting 'cosine_similarity')
    k_neighbors : int, optional
        number of k nearest neighbors to find between datasets
    consider_clusters : bool, optional
        flag indicating whether to only identify similar cells between datasets within the same annotated cluster
    save_chunks : bool, optional
        flag indicating whether to write intermediate results to the disk
    chunk_size : int, optional
        number of cells per chunk to impute at a time
    """

    x, y, clusters, save_location = read_inputs(
        matrix_to_impute_for=matrix_to_impute_for,
        matrix_to_impute_from=matrix_to_impute_from,
        cell_identities=cell_identities,
        cell_name_column_idx=cell_name_column_idx,
        outdir=outdir
    )

    similarity_output = similarity_matrix(
        x,
        y,
        clusters = clusters,
        cell_name_column_idx = cell_name_column_idx,
        metric = metric,
        k_neighbors = k_neighbors,
        consider_clusters = consider_clusters,
        outdir = save_location
    )

    imputed_results = impute_expression(
        x,
        y,
        similarity_matrix=similarity_output,
        chunk_size = chunk_size,
        outdir = save_location
    )

    final_df = merge_imputation_chunks(
        imputations_dir = save_location,
        outdir = save_location,
        save_chunks = save_chunks
    )

    validate_results(
        x = x,
        imputed_mtx = final_df,
        outdir = save_location
    )