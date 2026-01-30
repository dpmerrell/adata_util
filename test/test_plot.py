"""Tests for the plot command."""

import os
import subprocess

import anndata as ad
import numpy as np
import pandas as pd


def test_plot_embedding_pca(temp_dir):
    """Test plotting PCA embedding."""
    np.random.seed(42)
    n_obs, n_vars = 50, 20
    X = np.random.randn(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame(
        {"cluster": np.random.choice(["A", "B", "C"], n_obs)},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    adata = ad.AnnData(X=X, obs=obs)
    adata.obsm["X_pca"] = np.random.randn(n_obs, 10).astype(np.float32)

    input_path = os.path.join(temp_dir, "with_pca.h5ad")
    adata.write_h5ad(input_path)

    output_path = os.path.join(temp_dir, "pca_plot.png")

    result = subprocess.run(
        ["adata_util", "plot", "embedding", input_path, output_path, "--basis", "X_pca"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)
    assert os.path.getsize(output_path) > 0


def test_plot_embedding_with_color(temp_dir):
    """Test plotting embedding with color kwarg."""
    np.random.seed(42)
    n_obs, n_vars = 50, 20
    X = np.random.randn(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame(
        {"cluster": np.random.choice(["A", "B", "C"], n_obs)},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    adata = ad.AnnData(X=X, obs=obs)
    adata.obsm["X_pca"] = np.random.randn(n_obs, 10).astype(np.float32)

    input_path = os.path.join(temp_dir, "with_pca.h5ad")
    adata.write_h5ad(input_path)

    output_path = os.path.join(temp_dir, "pca_colored.png")

    result = subprocess.run(
        [
            "adata_util",
            "plot",
            "embedding",
            input_path,
            output_path,
            "--basis",
            "X_pca",
            "--color",
            "cluster",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "color" in result.stdout
    assert os.path.exists(output_path)


def test_plot_embedding_pdf(temp_dir):
    """Test plotting embedding to PDF."""
    np.random.seed(42)
    n_obs, n_vars = 50, 20
    X = np.random.randn(n_obs, n_vars).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obsm["X_pca"] = np.random.randn(n_obs, 10).astype(np.float32)

    input_path = os.path.join(temp_dir, "with_pca.h5ad")
    adata.write_h5ad(input_path)

    output_path = os.path.join(temp_dir, "pca_plot.pdf")

    result = subprocess.run(
        ["adata_util", "plot", "embedding", input_path, output_path, "--basis", "X_pca"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)
    assert os.path.getsize(output_path) > 0


def test_plot_embedding_with_title(temp_dir):
    """Test plotting embedding with title kwarg."""
    np.random.seed(42)
    n_obs = 50
    X = np.random.randn(n_obs, 20).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obsm["X_pca"] = np.random.randn(n_obs, 10).astype(np.float32)

    input_path = os.path.join(temp_dir, "with_pca.h5ad")
    adata.write_h5ad(input_path)

    output_path = os.path.join(temp_dir, "pca_titled.png")

    result = subprocess.run(
        [
            "adata_util",
            "plot",
            "embedding",
            input_path,
            output_path,
            "--basis",
            "X_pca",
            "--title",
            "My PCA Plot",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)


def test_plot_violin(temp_dir):
    """Test plotting violin plot."""
    np.random.seed(42)
    n_obs, n_vars = 50, 20
    X = np.random.randn(n_obs, n_vars).astype(np.float32)
    obs = pd.DataFrame(
        {"cluster": np.random.choice(["A", "B", "C"], n_obs)},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    input_path = os.path.join(temp_dir, "for_violin.h5ad")
    adata.write_h5ad(input_path)

    output_path = os.path.join(temp_dir, "violin.png")

    result = subprocess.run(
        [
            "adata_util",
            "plot",
            "violin",
            input_path,
            output_path,
            "--keys",
            "gene_0",
            "gene_1",
            "--groupby",
            "cluster",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)
    assert os.path.getsize(output_path) > 0
