import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
from pathlib import Path
from scimpute.utils import intersection, checkformat
from scimpute.io import read_matrix, read_matrix_basic


def validate_results(
    x,
    imputed_mtx,
    outdir
):
    """
    Quality check of the gene expression imputation results by comparing experimentally measured values to imputed values.

    Parameters
    ----------
    x : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to the gene expression matrix to perform gene expression imputation for (e.g. from a spatial experiment); genes as columns, cells/spots as rows
    imputed_mtx : str or pathlib.Path or pandas.DataFrame, required
        pandas.DataFrame or path to the imputed gene expression matrix
    outdir : str or pathlib.Path, required
        location of the output directory
    """

    xformat = checkformat(x)
    imputed_mtxformat = checkformat(imputed_mtx)

    if xformat == "path":
        x = read_matrix(x)
    if imputed_mtxformat == "path":
        imputed_mtx = read_matrix_basic(imputed_mtx)



    x_true = imputed_mtx[imputed_mtx.columns.intersection(x.columns)]
    x_sorted = x[x_true.columns]


    qualities = []
    for i in range(len(x_true.index)):
        X = x_true.iloc[i].values
        Y = x_sorted.iloc[i].values
        X1 = X.reshape(1, -1)
        Y1 = Y.reshape(1, -1)

        cs = cosine_similarity(X1, Y1)
        qualities.append(cs[0,0])    
    plt.hist(qualities)


    fig, ax = plt.subplots(figsize=(15, 15))
    plt.rcParams["font.size"] = 18
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.hist(x=qualities, bins=40, range=(0,1))
    plt.xlabel("cosine similarity")
    plt.ylabel("frequency")
    plt.savefig(os.path.join(outdir, "cosine_similarities.png"), bbox_inches="tight")


    scores = np.array(qualities)
    mean = np.mean(scores)
    median = np.median(scores)
    lower_quantile = np.percentile(scores, 25)
    upper_quantile = np.percentile(scores, 75)
    stdev = np.std(scores)

    f = open(os.path.join(outdir, "stats.txt"), "w+")
    f.write(f"Mean: {mean}\nMedian: {median}\nLower quantile: {lower_quantile}\nUpper quantile: {upper_quantile}\nStandard deviation: {stdev}")
    f.close()

    print_str = " ".join(map(str, qualities))
    with open(os.path.join(outdir, "scores.txt"), "w") as f:
        f.write(print_str)

    print("Mean:", mean)
    print("Median:", median)
    print("Lower quantile:", lower_quantile)
    print("Upper quantile:", upper_quantile)
    print("Standard deviation:", stdev)
