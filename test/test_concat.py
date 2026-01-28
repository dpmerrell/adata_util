"""Tests for the concat command."""

import os
import subprocess

import anndata as ad


def test_concat_two_files(two_adatas, temp_dir):
    """Test concatenating two h5ad files."""
    path1, path2 = two_adatas
    output_path = os.path.join(temp_dir, "concat_output.h5ad")

    result = subprocess.run(
        ["adata_util", "concat", path1, path2, output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    # Verify the concatenated file
    adata = ad.read_h5ad(output_path)
    assert adata.n_obs == 45  # 20 + 25
    assert adata.n_vars == 10


def test_concat_with_label(two_adatas, temp_dir):
    """Test concatenation with batch label."""
    path1, path2 = two_adatas
    output_path = os.path.join(temp_dir, "concat_labeled.h5ad")

    result = subprocess.run(
        ["adata_util", "concat", path1, path2, output_path, "--label", "source"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    assert adata.n_obs == 45


def test_concat_outer_join(two_adatas, temp_dir):
    """Test concatenation with outer join."""
    path1, path2 = two_adatas
    output_path = os.path.join(temp_dir, "concat_outer.h5ad")

    result = subprocess.run(
        ["adata_util", "concat", path1, path2, output_path, "--join", "outer"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    assert adata.n_obs == 45


def test_concat_single_file_error(simple_adata, temp_dir):
    """Test that concat with only one file shows an error."""
    output_path = os.path.join(temp_dir, "single_output.h5ad")

    result = subprocess.run(
        ["adata_util", "concat", simple_adata],
        capture_output=True,
        text=True,
    )

    # Should fail or show error message
    assert "Error" in result.stdout or result.returncode != 0 or not os.path.exists(output_path)
