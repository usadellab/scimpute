from scimpute.io import read_matrix
import pandas as pd
from pathlib import Path

def intersection(lst1, lst2):
    """
    Finds intersection of items between two lists.

    Parameters
    ----------
    lst1 : list
        First list containing items.
    lst2 : list
        Second list containing items.
    """

    return list(set(lst1) & set(lst2))

def checkformat(x):
    """
    Checks if input is string, pathlib.Path or pandas.DataFrame and handles file reading accordingly.

    Parameters
    ----------
    x : str or pathlib.Path or pandas.DataFrame
        the input object that should be checked for data type
    """

    if isinstance(x, (Path, str)):
        return "path"
    elif isinstance(x, pd.DataFrame):
        return "dataframe"
    else:
        raise ValueError(f"The input format for {x} is invalid. {x} must be of format str or pathlib.Path or pandas.DataFrame.")

def download_demo():
    """

    Returns
    -------
    matrix_to_impute_for : pandas.DataFrame
        pandas.DataFrame demo data from https://doi.org/10.60534/zfdth-2g147 containing the gene expression matrix to perform gene expression imputation for (e.g. from a spatial experiment); genes as columnr, cells/spots as rows.
    matrix_to_impute_from : pandas.DataFrame
        pandas.DataFrame demo data from https://doi.org/10.60534/zfdth-2g147 containing the comprehensive gene expression matrix to perform gene expression imputation from (e.g. from a single cell/nucleus RNA sequencing experiment); genes as columns, cells/spots as rows.
    cell_identities : pandas.DataFrame
        pandas.DataFrame demo data from https://doi.org/10.60534/zfdth-2g147 containing the identities for all cells from both datasets
    cell_name_column_idx : int
        index of the column containing the cell names (not the cell type clusters)
    """
    
    matrix_to_impute_for = pd.read_csv("https://git.nfdi4plants.org/usadellab/Barvista_ARC/-/raw/main/runs/r_CPM_matrices/D2-1_CPM_normalized_IM.tsv", sep="\t", index_col=0).transpose()
    matrix_to_impute_from = pd.read_csv("https://git.nfdi4plants.org/usadellab/Barvista_ARC/-/raw/main/runs/r_CPM_matrices/CPM_normalized_vSAM.tsv", sep="\t", index_col=0).transpose()
    cell_name_column_idx = 1
    cell_identities = pd.read_csv("https://git.nfdi4plants.org/usadellab/Barvista_ARC/-/raw/main/runs/r_clustering/resolve_D2-1/D2-1_clusters_1.3.tsv", sep="\t", index_col=cell_name_column_idx)
    cell_identities.index = cell_identities.index.str.replace("edgar", "vSAM")
    
    return (matrix_to_impute_for, matrix_to_impute_from, cell_identities, cell_name_column_idx)