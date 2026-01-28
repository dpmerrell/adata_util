"""Tests for the join command."""

import os
import subprocess

import anndata as ad
import pandas as pd


def test_join_obs_csv(simple_adata, metadata_csv, temp_dir):
    """Test joining CSV metadata to obs."""
    output_path = os.path.join(temp_dir, "joined_obs.h5ad")

    result = subprocess.run(
        ["adata_util", "join", simple_adata, metadata_csv, output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    adata = ad.read_h5ad(output_path)
    # Original columns plus new columns (score, label)
    assert "score" in adata.obs.columns
    assert "label" in adata.obs.columns
    # Original columns preserved
    assert "cell_type" in adata.obs.columns
    assert "batch" in adata.obs.columns


def test_join_var_tsv(simple_adata, metadata_tsv, temp_dir):
    """Test joining TSV metadata to var."""
    output_path = os.path.join(temp_dir, "joined_var.h5ad")

    result = subprocess.run(
        ["adata_util", "join", simple_adata, metadata_tsv, output_path, "--on", "var"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    adata = ad.read_h5ad(output_path)
    assert "pathway" in adata.var.columns
    # Original columns preserved
    assert "gene_name" in adata.var.columns


def test_join_with_right_on(simple_adata, temp_dir):
    """Test join with explicit right-on column."""
    # Create a CSV with a different column name for the key
    df = pd.DataFrame(
        {
            "some_score": [0.1, 0.2, 0.3],
            "the_cell_id": ["cell_0", "cell_1", "cell_2"],
        }
    )
    csv_path = os.path.join(temp_dir, "alt_metadata.csv")
    df.to_csv(csv_path, index=False)

    output_path = os.path.join(temp_dir, "joined_right_on.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "join",
            simple_adata,
            csv_path,
            output_path,
            "--right-on",
            "the_cell_id",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    assert "some_score" in adata.obs.columns


def test_join_left_partial_match(simple_adata, temp_dir):
    """Test left join with partial matching (NaN for non-matches)."""
    # Create a CSV with only a few matching cell IDs
    df = pd.DataFrame(
        {
            "cell_id": ["cell_0", "cell_1", "cell_999"],  # cell_999 doesn't exist
            "value": [1.0, 2.0, 3.0],
        }
    )
    csv_path = os.path.join(temp_dir, "partial_metadata.csv")
    df.to_csv(csv_path, index=False)

    output_path = os.path.join(temp_dir, "joined_left.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "join",
            simple_adata,
            csv_path,
            output_path,
            "--how",
            "left",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    # Left join preserves all original rows
    assert adata.n_obs == 50
    assert "value" in adata.obs.columns
    # Only cell_0 and cell_1 should have values, rest should be NaN
    assert adata.obs.loc["cell_0", "value"] == 1.0
    assert adata.obs.loc["cell_1", "value"] == 2.0
    assert pd.isna(adata.obs.loc["cell_2", "value"])


def test_join_inner(simple_adata, temp_dir):
    """Test inner join subsets AnnData to matching keys only."""
    # Create a CSV with only a few matching cell IDs
    df = pd.DataFrame(
        {
            "cell_id": ["cell_0", "cell_1", "cell_999"],  # cell_999 doesn't exist
            "value": [1.0, 2.0, 3.0],
        }
    )
    csv_path = os.path.join(temp_dir, "partial_metadata.csv")
    df.to_csv(csv_path, index=False)

    output_path = os.path.join(temp_dir, "joined_inner.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "join",
            simple_adata,
            csv_path,
            output_path,
            "--how",
            "inner",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0

    adata = ad.read_h5ad(output_path)
    # Inner join keeps only the intersection of keys
    assert adata.n_obs == 2
    assert set(adata.obs_names) == {"cell_0", "cell_1"}
    assert "value" in adata.obs.columns
    assert adata.obs.loc["cell_0", "value"] == 1.0
    assert adata.obs.loc["cell_1", "value"] == 2.0
