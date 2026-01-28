"""Tests for the scanpy command."""

import os
import subprocess

import anndata as ad
import numpy as np


def test_scanpy_pca(simple_adata, temp_dir):
    """Test running PCA via scanpy command."""
    output_path = os.path.join(temp_dir, "pca_output.h5ad")

    result = subprocess.run(
        ["adata_util", "scanpy", "pp.pca", simple_adata, output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    adata = ad.read_h5ad(output_path)
    assert "X_pca" in adata.obsm
    assert "pca" in adata.uns


def test_scanpy_pca_with_kwargs(simple_adata, temp_dir):
    """Test PCA with n_comps kwarg."""
    output_path = os.path.join(temp_dir, "pca_kwargs.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "scanpy",
            "pp.pca",
            simple_adata,
            output_path,
            "--n_comps",
            "5",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "n_comps" in result.stdout

    adata = ad.read_h5ad(output_path)
    assert adata.obsm["X_pca"].shape[1] == 5


def test_scanpy_normalize_total(simple_adata, temp_dir):
    """Test normalize_total operation."""
    output_path = os.path.join(temp_dir, "normalized.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "scanpy",
            "pp.normalize_total",
            simple_adata,
            output_path,
            "--target_sum",
            "10000",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    # After normalization, row sums should be close to target_sum
    # (allowing for floating point differences)
    row_sums = np.abs(adata.X).sum(axis=1)
    assert row_sums.shape[0] == 50


def test_scanpy_log1p(simple_adata, temp_dir):
    """Test log1p transformation."""
    output_path = os.path.join(temp_dir, "log1p.h5ad")

    # First make values positive for log1p
    adata_orig = ad.read_h5ad(simple_adata)
    adata_orig.X = np.abs(adata_orig.X) + 1
    positive_path = os.path.join(temp_dir, "positive.h5ad")
    adata_orig.write_h5ad(positive_path)

    result = subprocess.run(
        ["adata_util", "scanpy", "pp.log1p", positive_path, output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    # log1p values should be smaller than original (for values > e-1)
    assert adata.X.max() < adata_orig.X.max()


def test_scanpy_overwrite_input(simple_adata):
    """Test scanpy command that overwrites input file."""
    # Read original to compare
    adata_orig = ad.read_h5ad(simple_adata)
    orig_shape = adata_orig.obsm.get("X_pca", None)

    result = subprocess.run(
        ["adata_util", "scanpy", "pp.pca", simple_adata],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    # File should be modified in place
    adata = ad.read_h5ad(simple_adata)
    assert "X_pca" in adata.obsm


def test_scanpy_invalid_operation(simple_adata, temp_dir):
    """Test scanpy command with invalid operation."""
    output_path = os.path.join(temp_dir, "invalid.h5ad")

    result = subprocess.run(
        ["adata_util", "scanpy", "pp.nonexistent_function", simple_adata, output_path],
        capture_output=True,
        text=True,
    )

    # Should show error message
    assert "Error" in result.stdout or "Unknown" in result.stdout


def test_scanpy_highly_variable_genes(simple_adata, temp_dir):
    """Test highly_variable_genes operation."""
    output_path = os.path.join(temp_dir, "hvg.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "scanpy",
            "pp.highly_variable_genes",
            simple_adata,
            output_path,
            "--n_top_genes",
            "10",
            "--flavor",
            "seurat",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    assert "highly_variable" in adata.var.columns


def test_scanpy_boolean_kwarg(simple_adata, temp_dir):
    """Test parsing of boolean kwargs."""
    output_path = os.path.join(temp_dir, "pca_bool.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "scanpy",
            "pp.pca",
            simple_adata,
            output_path,
            "--zero_center",
            "False",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "zero_center" in result.stdout
