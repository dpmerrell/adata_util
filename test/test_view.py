"""Tests for the view command."""

import subprocess


def test_view_simple(simple_adata):
    """Test view command on a simple AnnData."""
    result = subprocess.run(
        ["adata_util", "view", simple_adata],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    output = result.stdout

    # Check basic structure is reported
    assert "50 observations x 20 variables" in output
    assert "obs (observation/cell metadata): 2 columns" in output
    assert "cell_type" in output
    assert "batch" in output
    assert "var (variable/gene metadata): 1 columns" in output
    assert "gene_name" in output


def test_view_with_layers(adata_with_layers):
    """Test view command on AnnData with layers, obsm, uns."""
    result = subprocess.run(
        ["adata_util", "view", adata_with_layers],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    output = result.stdout

    # Check shape
    assert "30 observations x 15 variables" in output

    # Check layers reported
    assert "layers: 1 entries" in output
    assert "raw_counts" in output

    # Check obsm reported
    assert "obsm (observation matrices): 1 entries" in output
    assert "X_pca" in output

    # Check uns reported
    assert "uns (unstructured): 1 entries" in output
    assert "metadata" in output


def test_view_nonexistent_file(temp_dir):
    """Test view command with nonexistent file."""
    result = subprocess.run(
        ["adata_util", "view", f"{temp_dir}/nonexistent.h5ad"],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
