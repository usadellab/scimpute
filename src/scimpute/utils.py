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