# `adata_util`

A Python package and command line interface for operating on AnnData objects (h5ad files).

A thin wrapper around [`scanpy`](https://scanpy.readthedocs.io/) and [`anndata`](https://anndata.readthedocs.io/).

## Installation

**Requirements:** Python >= 3.13

Clone this repository and install via pip:

```bash
git clone git@github.com:dpmerrell/adata_util.git
cd adata_util
pip install -e .
```

### Tab Completion (Optional)

Enable tab completion for subcommands and arguments:

```bash
# For bash (add to ~/.bashrc)
eval "$(register-python-argcomplete adata_util)"

# For zsh (add to ~/.zshrc)
eval "$(register-python-argcomplete adata_util)"
```

## Commands

| Command | Description |
|---------|-------------|
| `view` | Print a summary of an h5ad file's contents |
| `concat` | Concatenate multiple h5ad files |
| `join` | Join tabular data to obs or var |
| `split` | Extract obs or var table to CSV/TSV |
| `scanpy` | Run any scanpy operation |
| `plot` | Run any scanpy plotting function |
| `score_genes` | Score gene lists using scanpy.tl.score_genes |

## Command Reference

### `view`

Print an informative summary of an h5ad file, including shape, metadata columns, embeddings, and sample values from X and layers.

```bash
adata_util view <input.h5ad>
```

**Example output:**
```
AnnData object: my_data.h5ad
  Shape: 1000 observations x 2000 variables

X (main data matrix):
    Storage: sparse, dtype: float32
    Sample nonzero values: 1, 3, 2, 5, 1, 4, 2, 1, 3, 2

obs (observation/cell metadata): 3 columns
    cell_type: category
    batch: category
    n_counts: float64
...
```

### `concat`

Concatenate AnnData objects from multiple h5ad files.

```bash
adata_util concat <input1.h5ad> <input2.h5ad> [...] <output.h5ad> [options]
```

**Options:**
- `--axis {0,1}`: Concatenate along obs (0) or var (1). Default: 0
- `--join {inner,outer}`: How to handle the other axis. Default: inner
- `--label <name>`: Add a column with batch labels

**Example:**
```bash
# Concatenate three samples
adata_util concat sample1.h5ad sample2.h5ad sample3.h5ad merged.h5ad

# Concatenate with batch labels
adata_util concat sample1.h5ad sample2.h5ad merged.h5ad --label batch
```

### `join`

Join obs- or var-level data from a tabular file (CSV/TSV) to an AnnData.

```bash
adata_util join <input.h5ad> <table.csv> <output.h5ad> [options]
```

**Options:**
- `--on {obs,var}`: Join to obs or var. Default: obs
- `--left-on <col>`: Column in AnnData to join on. Default: index
- `--right-on <col>`: Column in table to join on. Default: first column
- `--how {left,inner}`: Type of join. Default: left

**Example:**
```bash
# Add cell annotations from a CSV
adata_util join data.h5ad annotations.csv annotated.h5ad

# Inner join (keep only matching cells)
adata_util join data.h5ad metadata.tsv filtered.h5ad --how inner
```

### `split`

Extract the obs or var table from an AnnData to a tabular file.

```bash
adata_util split <input.h5ad> {obs,var} <output.csv> [options]
```

**Options:**
- `--kept-cols <col1> [col2] ...`: Columns to keep. Default: all columns

**Example:**
```bash
# Export all cell metadata
adata_util split data.h5ad obs cell_metadata.csv

# Export specific columns
adata_util split data.h5ad obs clusters.tsv --kept-cols cluster cell_type
```

### `scanpy`

Run a scanpy operation on an AnnData. Supports any scanpy function that operates on a single AnnData and modifies it in place.

```bash
adata_util scanpy <operation> <input.h5ad> [output.h5ad] [--kwargs]
```

**Arguments:**
- `operation`: Scanpy function (e.g., `pp.normalize_total`, `tl.pca`, `tl.umap`)
- `output`: Output file. If omitted, overwrites input.

**Example:**
```bash
# Standard preprocessing pipeline
adata_util scanpy pp.filter_cells raw.h5ad filtered.h5ad --min_genes 200
adata_util scanpy pp.filter_genes filtered.h5ad filtered.h5ad --min_cells 3
adata_util scanpy pp.normalize_total filtered.h5ad normalized.h5ad --target_sum 10000
adata_util scanpy pp.log1p normalized.h5ad normalized.h5ad
adata_util scanpy pp.highly_variable_genes normalized.h5ad normalized.h5ad --n_top_genes 2000
adata_util scanpy pp.pca normalized.h5ad processed.h5ad --n_comps 50
adata_util scanpy pp.neighbors processed.h5ad processed.h5ad
adata_util scanpy tl.umap processed.h5ad processed.h5ad
adata_util scanpy tl.leiden processed.h5ad processed.h5ad --resolution 0.5
```

### `plot`

Run a scanpy plotting function on an h5ad file. Supports any `scanpy.pl` function.

```bash
adata_util plot <function> <input.h5ad> <output.png> [--kwargs]
```

**Arguments:**
- `function`: Plotting function (e.g., `embedding`, `violin`, `heatmap`)
- `output`: Image file (.png, .pdf, .svg)

**Example:**
```bash
# Basic UMAP plot
adata_util plot embedding data.h5ad umap.png --basis X_umap

# Colored by cluster and gene expression
adata_util plot embedding data.h5ad umap_colored.png --basis X_umap --color leiden cell_type CD4

# Customized PCA plot
adata_util plot embedding data.h5ad pca.pdf --basis X_pca --color batch --palette Set2 --size 20

# Violin plot
adata_util plot violin data.h5ad violin.png --keys gene1 gene2 --groupby cluster

# Heatmap
adata_util plot heatmap data.h5ad heatmap.png --var_names gene1 gene2 gene3 --groupby cluster
```

### `score_genes`

Score gene lists using scanpy.tl.score_genes. Useful for scoring pathway activity or cell type signatures.

```bash
adata_util score_genes <input.h5ad> <genelists.yaml> <output.h5ad> [--kwargs]
```

**Arguments:**
- `genelists`: YAML or JSON file with gene lists (dict of name → gene list)

**Example gene list file (genelists.yaml):**
```yaml
cell_cycle:
  - CDK1
  - CCNB1
  - CCNA2
apoptosis:
  - BAX
  - BCL2
  - CASP3
```

**Example:**
```bash
adata_util score_genes data.h5ad signatures.yaml scored.h5ad

# With custom control size
adata_util score_genes data.h5ad signatures.json scored.h5ad --ctrl_size 100
```

This adds `cell_cycle_score` and `apoptosis_score` columns to obs.

## Workflow Example

A complete single-cell RNA-seq analysis workflow:

```bash
# 1. View raw data
adata_util view raw_counts.h5ad

# 2. Basic QC and preprocessing
adata_util scanpy pp.filter_cells raw_counts.h5ad qc.h5ad --min_genes 200
adata_util scanpy pp.filter_genes qc.h5ad qc.h5ad --min_cells 3
adata_util scanpy pp.normalize_total qc.h5ad normalized.h5ad --target_sum 10000
adata_util scanpy pp.log1p normalized.h5ad normalized.h5ad

# 3. Feature selection and dimensionality reduction
adata_util scanpy pp.highly_variable_genes normalized.h5ad hvg.h5ad --n_top_genes 2000
adata_util scanpy pp.pca hvg.h5ad pca.h5ad --n_comps 50

# 4. Clustering
adata_util scanpy pp.neighbors pca.h5ad neighbors.h5ad
adata_util scanpy tl.umap neighbors.h5ad umap.h5ad
adata_util scanpy tl.leiden umap.h5ad clustered.h5ad --resolution 0.8

# 5. Score gene signatures
adata_util score_genes clustered.h5ad cell_signatures.yaml scored.h5ad

# 6. Add external annotations
adata_util join scored.h5ad cell_annotations.csv annotated.h5ad

# 7. Generate plots
adata_util plot embedding annotated.h5ad umap_clusters.png --basis X_umap --color leiden
adata_util plot embedding annotated.h5ad umap_celltypes.png --basis X_umap --color cell_type

# 8. Export metadata for downstream analysis
adata_util split annotated.h5ad obs cell_metadata.csv
```

## License

MIT
