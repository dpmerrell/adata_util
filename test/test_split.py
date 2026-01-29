"""Tests for the split command."""

import os
import subprocess

import pandas as pd


def test_split_obs_csv(simple_adata, temp_dir):
    """Test splitting obs table to CSV."""
    output_path = os.path.join(temp_dir, "obs.csv")

    result = subprocess.run(
        ["adata_util", "split", simple_adata, "obs", output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    df = pd.read_csv(output_path, index_col=0)
    assert "cell_type" in df.columns
    assert "batch" in df.columns
    assert len(df) == 50


def test_split_var_tsv(simple_adata, temp_dir):
    """Test splitting var table to TSV."""
    output_path = os.path.join(temp_dir, "var.tsv")

    result = subprocess.run(
        ["adata_util", "split", simple_adata, "var", output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    df = pd.read_csv(output_path, sep="\t", index_col=0)
    assert "gene_name" in df.columns
    assert len(df) == 20


def test_split_kept_cols(simple_adata, temp_dir):
    """Test splitting with specific columns kept."""
    output_path = os.path.join(temp_dir, "obs_subset.csv")

    result = subprocess.run(
        [
            "adata_util",
            "split",
            simple_adata,
            "obs",
            output_path,
            "--kept-cols",
            "cell_type",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    df = pd.read_csv(output_path, index_col=0)
    assert "cell_type" in df.columns
    assert "batch" not in df.columns
    assert len(df.columns) == 1


def test_split_kept_cols_multiple(simple_adata, temp_dir):
    """Test splitting with multiple columns kept."""
    output_path = os.path.join(temp_dir, "obs_multi.csv")

    result = subprocess.run(
        [
            "adata_util",
            "split",
            simple_adata,
            "obs",
            output_path,
            "--kept-cols",
            "cell_type",
            "batch",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    df = pd.read_csv(output_path, index_col=0)
    assert "cell_type" in df.columns
    assert "batch" in df.columns
    assert len(df.columns) == 2


def test_split_missing_column(simple_adata, temp_dir):
    """Test error when specifying nonexistent column."""
    output_path = os.path.join(temp_dir, "obs_error.csv")

    result = subprocess.run(
        [
            "adata_util",
            "split",
            simple_adata,
            "obs",
            output_path,
            "--kept-cols",
            "nonexistent_column",
        ],
        capture_output=True,
        text=True,
    )

    assert "Error" in result.stdout or "not found" in result.stdout


def test_split_index_preserved(simple_adata, temp_dir):
    """Test that index is always preserved."""
    output_path = os.path.join(temp_dir, "obs_index.csv")

    result = subprocess.run(
        [
            "adata_util",
            "split",
            simple_adata,
            "obs",
            output_path,
            "--kept-cols",
            "cell_type",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    df = pd.read_csv(output_path, index_col=0)
    # Index should be cell names
    assert "cell_0" in df.index
    assert "cell_49" in df.index
