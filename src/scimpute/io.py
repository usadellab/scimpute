import pandas as pd
from pathlib import Path
import os

def read_inputs(
    matrix_to_impute_for,
    matrix_to_impute_from,
    cell_identities,
    cell_name_column_idx,
    outdir
):
    """
    Reading all needed input files.

    Parameters
    ----------
    matrix_to_impute_for : str or pathlib.Path, required
        gene expression matrix to perform gene expression imputation for (e.g. from a spatial experiment)
    matrix_to_impute_from : str or Path, required
        comprehensive gene expression matrix to perform gene expression imputation from (e.g. from a single cell/nucleus RNA sequencing experiment)
    cell_identities : str or Path, required
        table containing identities for all cells from both datasets
    cell_name_column_idx : int, required
        index of the column containing the cell names (not the cell type clusters)
    outdir : str or pathlib.Path, required
        location of output directory
    """

    if os.path.isabs(outdir):
        save_location = outdir
    else:
        save_location = os.path.join(os.getcwd(), outdir)
    Path(save_location).mkdir(parents=True, exist_ok=True)

    x = read_matrix(matrix_to_impute_for)
    y = read_matrix(matrix_to_impute_from)

    clusters = read_cell_identities(cell_identities, index_column=cell_name_column_idx)

    return (x, y, clusters, Path(save_location))

def read_matrix(
    filepath,
    separator="\t"
):
    """
    Reads a gene expression matrix.

    Parameters
    ----------
    filepath : str or pathlib.Path, required
        filepath to the gene expression file
    separator : str or RegEx, optional
        field separator of the gene expression file
    """

    mtx = Path(filepath)
    if not mtx.exists():
        raise FileNotFoundError(f"Matrix file not found: {mtx}")
    mtx = pd.read_csv(filepath, sep=separator, index_col=0)
    mtx = mtx.transpose()
    
    return mtx

def read_cell_identities(
    filepath,
    index_column,
    separator="\t"
):
    """
    Reads a cell type identity file.

    Parameters
    ----------
    filepath : str or pathlib.Path, required
        filepath to the cell type identity file
    index_column : int, required
        index of the column which contains the cell names (not the cell type clusters)
    separator : str or RegEx, optional
        field separator of the gene expression file
    """

    clusters = Path(filepath)
    if not clusters.exists():
        raise FileNotFoundError(f"Cell type identity file not found: {clusters}")
    clusters = pd.read_csv(filepath, sep=separator, index_col=index_column)
    clusters.index = clusters.index.str.lstrip("_")
    clusters.rename(columns={
        clusters.columns[0]: "id"
    })

    return clusters

def read_matrix_basic(
    filepath,
    separator="\t"
):
    """
    Reads a table.

    Parameters
    ----------
    filepath : str or pathlib.Path, required
        filepath to the cell type identity file
    separator : str or RegEx, optional
        field separator of the gene expression file
    """

    mtx = pd.read_csv(filepath, sep=separator)
    return mtx