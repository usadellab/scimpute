import os
import pandas as pd
import shutil
from pathlib import Path

def merge_imputation_chunks(
    imputations_dir,
    outdir,
    save_chunks = False
):
    """
    Merges multiple imputed gene expression matrices from within one folder.

    Parameters
    ----------
    imputations_dir : str or pathlib.Path, required
        path to the directory containing the imputation files
    outdir : str or pathlib.Path, required
        output directory for the final imputed gene expression matrix
    save_chunks : bool, optional
        flag indicating whether to write keep unmerged results on the disk
    """

    Path(outdir).mkdir(parents=True, exist_ok=True)

    if imputations_dir != outdir:
        chunkdir = imputations_dir
    else:
        chunkdir = os.path.join(outdir, "imputations")

    files = [os.path.join(chunkdir, file) for file in os.listdir(chunkdir)]
    files = sorted(files)

    final_df = pd.DataFrame()
    for file in files:
        print(f"Merging file {file}")
        df = pd.read_csv(file, sep="\t", index_col=0)
        final_df = pd.concat([final_df, df])

    final_df.index = final_df["cell"]
    final_df.drop(columns=["cell"], inplace=True)
    final_df.to_csv(os.path.join(outdir, "imputation.tsv"), sep="\t")

    if not save_chunks:
        shutil.rmtree(chunkdir)
    
    return final_df