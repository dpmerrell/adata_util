# `adata_util`

A python package and command line interface for operating on AnnData objects (i.e., h5ad files).

A thin wrapper around [`scanpy`](https://scanpy.readthedocs.io/).

# Installation

Download or clone this repository, and then install the package via pip: 
```bash
git clone git@github.com:dpmerrell/adata_util.git
cd adata_util
pip install -e .
```

# Basic usage

The `adata_util` package comes with the `adata_util` terminal command, which possesses a large collection of subcommands:

* `view`: prints an informative summary of the contents of an h5ad file
* `concat`: concatenates the AnnData objects in one or more h5ad files
* `join`: join obs- or var- level data from a tabular file (e.g., tsv) to an AnnData in an h5ad file.
* `scanpy`: specify a `scanpy` operation to run on the AnnData in an h5ad file. (Restricted to `scanpy` operations that (i) operate on a single AnnData object and (ii) append results to that AnnData).
* `plot_embedding`: plot an embedding (e.g., PCA, UMAP) from an h5ad file to an image file. A thin wrapper around `scanpy.pl.embedding`.

Most of these commands should output a modified version of the input h5ad(s).
That is, the user should specify an output h5ad that, in general, will differ from the input h5ad(s).
In other words, it shouldn't modify the input h5ad(s) by default.

## Examples:

```bash
$ adata_util view my_adata.h5ad
# <outputs a summary of the AnnData to stdout>

# Concatenate three AnnDatas and output to a new h5ad (my_adata_concat.h5ad)
$ adata_util concat my_adata_1.h5ad my_adata_2.h5ad my_adata_3.h5ad my_adata_concat.h5ad

# Apply some basic scanpy operations to the h5ad
$ adata_util scanpy pp.highly_variable_genes my_adata_concat.h5ad temp.h5ad
$ adata_util scanpy pp.normalize_total temp.h5ad temp.h5ad --target_sum 10000
$ adata_util scanpy pp.log1p temp.h5ad temp.h5ad
$ adata_util scanpy pp.pca temp.h5ad

# Plot the PCA embedding, colored by a metadata column
$ adata_util plot_embedding temp.h5ad X_pca pca_plot.png --color cluster
```

