"""Pytest fixtures for adata_util tests."""

import os
import tempfile
import shutil

import anndata as ad
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    dirpath = tempfile.mkdtemp()
    yield dirpath
    shutil.rmtree(dirpath)


@pytest.fixture
def simple_adata(temp_dir):
    """Create a simple AnnData object and save it to h5ad."""
    np.random.seed(42)
    n_obs, n_vars = 50, 20

    X = np.random.randn(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cell_type": np.random.choice(["TypeA", "TypeB", "TypeC"], n_obs),
            "batch": np.random.choice(["batch1", "batch2"], n_obs),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"gene_name": [f"gene_{i}" for i in range(n_vars)]},
        index=[f"gene_{i}" for i in range(n_vars)],
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    filepath = os.path.join(temp_dir, "simple.h5ad")
    adata.write_h5ad(filepath)

    return filepath


@pytest.fixture
def adata_with_layers(temp_dir):
    """Create an AnnData with layers, obsm, and uns."""
    np.random.seed(123)
    n_obs, n_vars = 30, 15

    X = np.random.randn(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame(
        {"cluster": np.random.choice(["A", "B"], n_obs)},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"highly_variable": np.random.choice([True, False], n_vars)},
        index=[f"gene_{i}" for i in range(n_vars)],
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["raw_counts"] = np.random.randint(0, 100, (n_obs, n_vars))
    adata.obsm["X_pca"] = np.random.randn(n_obs, 5).astype(np.float32)
    adata.uns["metadata"] = {"version": "1.0"}

    filepath = os.path.join(temp_dir, "with_layers.h5ad")
    adata.write_h5ad(filepath)

    return filepath


@pytest.fixture
def two_adatas(temp_dir):
    """Create two AnnData objects for concatenation tests."""
    np.random.seed(456)

    # First adata
    X1 = np.random.randn(20, 10).astype(np.float32)
    obs1 = pd.DataFrame(
        {"sample": ["sample1"] * 20},
        index=[f"cell_a_{i}" for i in range(20)],
    )
    var = pd.DataFrame(
        {"gene_name": [f"gene_{i}" for i in range(10)]},
        index=[f"gene_{i}" for i in range(10)],
    )
    adata1 = ad.AnnData(X=X1, obs=obs1, var=var)
    path1 = os.path.join(temp_dir, "adata1.h5ad")
    adata1.write_h5ad(path1)

    # Second adata
    X2 = np.random.randn(25, 10).astype(np.float32)
    obs2 = pd.DataFrame(
        {"sample": ["sample2"] * 25},
        index=[f"cell_b_{i}" for i in range(25)],
    )
    adata2 = ad.AnnData(X=X2, obs=obs2, var=var.copy())
    path2 = os.path.join(temp_dir, "adata2.h5ad")
    adata2.write_h5ad(path2)

    return path1, path2


@pytest.fixture
def metadata_csv(temp_dir):
    """Create a CSV file with metadata to join."""
    df = pd.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(50)],
            "score": np.random.rand(50),
            "label": np.random.choice(["pos", "neg"], 50),
        }
    )
    filepath = os.path.join(temp_dir, "metadata.csv")
    df.to_csv(filepath, index=False)
    return filepath


@pytest.fixture
def metadata_tsv(temp_dir):
    """Create a TSV file with metadata to join."""
    df = pd.DataFrame(
        {
            "gene_id": [f"gene_{i}" for i in range(20)],
            "pathway": np.random.choice(["pathwayA", "pathwayB", "pathwayC"], 20),
        }
    )
    filepath = os.path.join(temp_dir, "metadata.tsv")
    df.to_csv(filepath, sep="\t", index=False)
    return filepath
