# Imputation of single cell/nucleus gene expression data

**scimpute** is a Python module for running gene expression imputation for single cell/nucleus RNA sequencing data. It has been developed within the frame of a project investigating development of the barley (*Hordeum vulgare*) shoot meristem, and was first applied in the corresponding publication of [Demesa-Arevalo et al. 2026](https://doi.org/10.1038/s41477-025-02176-6)[^1]. For ease of use, the gene expression imputation method is now available as a python module.

- [Imputation of single cell/nucleus gene expression data](#imputation-of-single-cellnucleus-gene-expression-data)
  - [Installation](#installation)
    - [Dependencies](#dependencies)
    - [User installation](#user-installation)
    - [Installation verification with demo data](#installation-verification-with-demo-data)
  - [Quickstart](#quickstart)
    - [Input files](#input-files)
    - [Single-command execution](#single-command-execution)
      - [Parameters](#parameters)
        - [Required parameters](#required-parameters)
        - [Optional parameters:](#optional-parameters)
    - [Step-wise execution](#step-wise-execution)
  - [Output files](#output-files)
  - [Concept](#concept)
  - [References](#references)



## Installation

### Dependencies

**scimpute** was developed using the following dependencies:
- Python >= 3.11.3
- setuptools >= 67.6.1
- numpy >= 1.24.2
- matplotlib >= 3.7.1
- scikit-learn >= 1.2.2

Using package versions higher than indicated here might lead to incompatibilities as some software may update in a way that it is not backwards compatible. In case of errors, please install the exact versions indicated here.

### User installation

**scimpute** can be easily installed using `pip`. First, clone this repository. Then, from the repository root, install **scimpute** along with all its dependencies:
```shell
git clone git@github.com:usadellab/scimpute.git
cd scimpute
pip install -e .
```

### Installation verification with demo data

To verify the installation on your machine, a test run with public demo data can be performed. The data (and the provided result data) belongs to the publication of [Demesa-Arevalo et al. 2026](https://doi.org/10.1038/s41477-025-02176-6)[^1], and is available at the corresponding data object at at [Dörpholz et al. 2025](https://doi.org/10.60534/zfdth-2g147)[^2]. To perform a test run of the installation with the demo data, the following python code can be executed, with a replacement of the `<output_path>` placeholder:

```python
from scimpute import expression_imputation
from scimpute.utils import download_demo

mtx_to_impute_for, mtx_to_impute_from, cl_identities, cl_name_column_idx = download_demo()

expression_imputation(
    # required parameters
    matrix_to_impute_for=mtx_to_impute_for,
    matrix_to_impute_from=mtx_to_impute_from,
    cell_identities=cl_identities,
    cell_name_column_idx=cl_name_column_idx,
    outdir=<output_path>,

    # optional parameters
    chunk_size=300,
    k_neighbors=4,
    consider_clusters=False,
    save_chunks=False
)
```

The result will be a directory in your specified `outdir` location, containing the five files detailed in the [output files](#output-files) section. You can compare the files you obtained with the demo outputs in the [tests](./tests/) directory of this repo.

## Quickstart

### Input files

To run the gene expression imputation on your data, the following input files are required:

- a gene expression matrix of the dataset you want to impute the gene expression for (e.g. a spatial transcriptomics matrix); tab-separated, genes as rows, cells as columns:
  ||ds1_cell1|ds1_cell2|ds1_cell3|
  |-|-|-|-|
  |gene1|0|40|0|
  |gene2|298|319|90|
- a gene expression matrix of a reference dataset, which you want to use to impute the gene expression from (e.g. a single cell RNAseq matrix); tab-separated, genes as rows, cells as columns:
  ||ds2_cell1|ds2_cell2|ds3_cell3|
  |-|-|-|-|
  |gene1|46512|0|0|
  |gene2|139535|210526|200000|
- a tab-separated table containing all cell names (as they are used in the two gene expression matrices) and their cluster identity as numbers, not as annotated label; the identity column must be named "id":
  |id|ident|
  |-|-|
  |2|ds1_cell1|
  |7|ds1_cell2|
  |2|ds2_cell1|

Ideally, the two gene expression matrices should be normalized/transformed in a comparable manner, such as using both matrices with CPM (counts per million) values, so that the similarity determination of cells between datasets can run smoothly (see section [concept](#concept)).

### Single-command execution

The entire pipeline can be executed using a single command. This is suitable if both datasets are small, or if enough RAM is available for calculation. Otherwise, please refer to the step-wise execution.

To run the gene expression imputation, the following command can be used after importing the **scimpute** package:

```python
from scimpute import expression_imputation

expression_imputation(
    # required parameters
    matrix_to_impute_for=<your_first_datset>,
    matrix_to_impute_from=<your_second_dataset>,
    cell_identities=<your_cell_identity_table>,
    cell_name_column_idx=<column_index_of_the_id_column>,
    outdir=<output_path>,

    # optional parameters
    chunk_size=1000,
    k_neighbors=25,
    consider_clusters=True,
    save_chunks=False
)
```
Replace the `<placeholder>` information with paths to your own files.

#### Parameters

:exclamation: On Windows, make sure that all slashes in your file paths are masked using another slash: e.g. `C:\User\Downloads` must be re-formatted to `C:\\User\\Downloads`.

##### Required parameters

- `matrix_to_impute_for`: enter the path to the gene expression matrix (in quotation marks) you want to perform the gene expression for (like a matrix from a spatial transcriptomics experiment with only few genes)
- `matrix_to_impute_from`: enter the path to the gene expression matrix (in quotation marks) you want to use as a reference to impute the gene expression from (likely a matrix from a single cell/nucleus RNAseq experiment with many genes)
- `cell_identities`: enter the path to the table containing the cluster identities for each cell (in quotation marks)
- `cell_name_column_idx`: enter the index of the column in the `cell_identities` table which contains the cell labels (index starts at 0)
- `outdir`: provide a path where your files should be saved (in quotation marks); if an absolute path is given, a folder will be created in this location; if a relative path is given, a folder will be created in your current working directory

##### Optional parameters:

- `chunk_size`: the number of cells per chunk to impute at a time; this is useful if the dataset is very large or if little RAM is available; default is `1000`
- `k_neighbors`: number of k nearest neighbors to find between datasets; default is `25`
- `consider_clusters`: a flag indicating whether to only identify similar cells between datasets within only the same annotated cluster or within the entire dataset; setting this to `True` will be less resource-intensive, but requires a reliably clustering to correctly identify the most similar cells; default is `True`
- `save_chunks`: flag indicating whether to keep intermediate results on the disk even after file merging; setting this to `True` is recommended for troubleshooting only; default is `False`


### Step-wise execution

The entire pipeline can be also be executed step by step. This is suitable if at least one of the datasets is large, or if only little RAM is available for calculation.

To run the gene expression imputation, the following commands can be used after importing individual modules of the **scimpute** package:

```python
from scimpute.io import read_inputs
from scimpute.similarity import similarity_matrix
from scimpute.imputation import impute_expression
from scimpute.merge import merge_imputation_chunks
from scimpute.validation import validate_results


x, y, clusters, save_location = read_inputs(
    matrix_to_impute_for=<your_first_dataset>, 
    matrix_to_impute_from=<your_second_dataset>, 
    cell_identities=<your_cell_identity_table>,
    cell_name_column_idx=<column_index_of_the_id_column>, 
    outdir=<output_path>
    )

similarity_output = similarity_matrix(
    x=<your_first_dataset>, 
    y=<your_second_dataset>,
    clusters=<your_cell_identity_table>,
    cell_name_column_idx=<column_index_of_the_id_column>, 
    outdir=<output_path>, 
    k_neighbors=<number_of_neighbors_you_want>, 
    consider_clusters=<flag_whether_you_want_to_consider_cluster_identity>
    )

imputed_results = impute_expression(
    x=<your_first_dataset>, 
    y=<your_second_dataset>, 
    similarity_matrix=<your_similarity_output_table_from_previous_step>, 
    outdir=<output_path>, 
    chunk_size=<number_of_cells_to_impute_for_at_a_time>
    )

final_df = merge_imputation_chunks(
    imputations_dir=<directory_containing_the_imputation_files>, 
    outdir=<output_path>, 
    save_chunks=<flag_whether_to_keep_unmerged_imputation_files>
    )

validate_results(
    x=<your_first_dataset>, 
    imputed_mtx=<your_imputed_gene_expression_matrix>, 
    outdir=<output_path>
    )
```

Please consult the [parameters](#parameters) section for details on the relevant parameters that can be used. In the step-by-step workflow, the `<placeholders>` for the following parameters must be given in quotation marks:
- `matrix_to_impute_for`
- `matrix_to_impute_from`
- `cell_identities`
- `outdir`
- `x`
- `y`
- `clusters`
- `similarity_matrix`
- `imputations_dir`
- `imputed_mtx`

## Output files

The tool will generate five final output files:
- `imputation.tsv`: the most important output file, which contains the imputed gene expression values for the dataset you wanted to run the gene expression for; the table is tab-separated, containing float values with two decimal points; genes are listed as columns, cells are listed as rows
  |cell|gene1|gene2|
  |-|-|-|
  |ds1_cell1|61.26|4.89|
  |ds1_cell2|16.96|0.0|
- `sim_matrix.tsv`: the results from the identification of the next most similar cells in the reference dataset for each cell in the dataset that you wanted to run the imputation for; the table is tab-separated, containing the cluster identity for each cell in the dataset you wanted to run the imputation for (column `cluster`), the name of the cell in this dataset (column `x`), the most similar cells in the reference dataset (column `y`), the cosine similarity between the cells in `x` and `y` (column `distance`) which is shifted by +1 into the positive range
  |cluster|x|y|distance|
  |-|-|-|-|
  |2|ds1_cell1|ds2_cell567|1.877341|
  |2|ds1_cell1|ds2_cell478|1.691159|
- `scores.txt`: a text file containing all similarity scores from the validation; each number indicates the cosine similarity for one cell in the dataset you wanted to run the imputation between the experimentally measured gene expression values and the imputed gene expression values; the number of scores is identical with the number of cells in your query dataset
- `cosine_similarity.png`: a histogram of the `scores.txt`, summarizing how well the imputation worked; a good imputation is indicated by a sharp peak at *1*, meaning that for most cells the gene expression pattern could be reproduced almost perfectly; a bad imputation is indicated by a broad distribution of scores across the x-axis
- `stats.txt`: a text file containing some statistical values of the cosine similarities displayed in the histogram: mean, Q25, Q50, Q75, standard deviation; the higher the mean and Q50, the better the imputation worked (this will be indicated by a sharp peak on the right-hand side in the histogram)

If the parameter `save_chunks=True` was set, there will be an additional folder which contains the same data as the `imputation.tsv`, but split into subsets of cells equivalent to the `chunk_size` (by default each submatrix contains imputed gene expression values for 1000 cells).

## Concept

The concept for the gene expression imputation procedure is briefly outlined here. For a more comprehensive explanation, please refer to [Demesa-Arevalo et al. 2026](https://doi.org/10.1038/s41477-025-02176-6)[^1]. Briefly, one dataset, usually a spatial transcriptomics dataset which contains information about few genes, is used as a query dataset to run the gene expression imputation for. A more comprehensive dataset, such as a single cell RNAseq dataset which contains information about many genes, is used as a reference dataset to compute gene expression from. It is assumed that the datasets were integrated and clustered together beforehand in order to identify which cells belong to the same cluster. For each cell in the query dataset, the *25* most similar cells in the reference dataset are identified based on the expression of the genes measured in both datasets, using cosine similarity as distance metric. From these most similar neighbors, a weighted average is calculated, using the cosine similarity as weight. This way, the gene expression values are computed for each gene measured in the reference dataset, and transferred to each cell in the query dataset. For internal validation, the imputed values for the limited number of genes measured in the query dataset are compared to the experimentally measured values for each cell in the dataset. This generates another cosine similarity score, which indicates whether the gene expression pattern could successfully be reproduced computationally. A cosine similarity score of *0* indicates no similarity, and a score of *+1* indicates perfect similarity.


## References

[^1]: Demesa-Arevalo, E., Dörpholz, H., Vardanega, I. *et al*. Imputation integrates single-cell and spatial gene expression data to resolve transcriptional networks in barley shoot meristem development. *Nat. Plants* (2026). https://doi.org/10.1038/s41477-025-02176-6
[^2]: Dörpholz, H., Demesa-Arevalo, E., Usadel, B., & Simon, R. (2025). Imputation integrates single-cell and spatial gene expression data to resolve transcriptional networks in barley shoot meristem development [Data set]. *DataPLANT*. https://doi.org/10.60534/zfdth-2g147