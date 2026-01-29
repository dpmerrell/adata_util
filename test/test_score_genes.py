"""Tests for the score_genes command."""

import json
import os
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
import yaml


def test_score_genes_json(temp_dir):
    """Test scoring gene lists from JSON file."""
    # Create test AnnData with enough genes for control gene selection
    np.random.seed(42)
    n_obs, n_vars = 50, 100
    X = np.abs(np.random.randn(n_obs, n_vars).astype(np.float32))
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, var=var)
    input_path = os.path.join(temp_dir, "input.h5ad")
    adata.write_h5ad(input_path)

    # Create gene lists JSON
    genelists = {
        "list_a": ["gene_0", "gene_1", "gene_2"],
        "list_b": ["gene_5", "gene_6"],
    }
    genelists_path = os.path.join(temp_dir, "genelists.json")
    with open(genelists_path, "w") as f:
        json.dump(genelists, f)

    output_path = os.path.join(temp_dir, "output.h5ad")

    result = subprocess.run(
        ["adata_util", "score_genes", input_path, genelists_path, output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    adata_out = ad.read_h5ad(output_path)
    assert "list_a_score" in adata_out.obs.columns
    assert "list_b_score" in adata_out.obs.columns


def test_score_genes_yaml(temp_dir):
    """Test scoring gene lists from YAML file."""
    # Create test AnnData with enough genes for control gene selection
    np.random.seed(42)
    n_obs, n_vars = 50, 100
    X = np.abs(np.random.randn(n_obs, n_vars).astype(np.float32))
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, var=var)
    input_path = os.path.join(temp_dir, "input.h5ad")
    adata.write_h5ad(input_path)

    # Create gene lists YAML
    genelists = {
        "pathway_x": ["gene_0", "gene_1"],
        "pathway_y": ["gene_10", "gene_11", "gene_12"],
    }
    genelists_path = os.path.join(temp_dir, "genelists.yaml")
    with open(genelists_path, "w") as f:
        yaml.dump(genelists, f)

    output_path = os.path.join(temp_dir, "output.h5ad")

    result = subprocess.run(
        ["adata_util", "score_genes", input_path, genelists_path, output_path],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert os.path.exists(output_path)

    adata_out = ad.read_h5ad(output_path)
    assert "pathway_x_score" in adata_out.obs.columns
    assert "pathway_y_score" in adata_out.obs.columns


def test_score_genes_with_kwargs(temp_dir):
    """Test scoring with additional kwargs."""
    # Create test AnnData with enough genes for control gene selection
    np.random.seed(42)
    n_obs, n_vars = 50, 100
    X = np.abs(np.random.randn(n_obs, n_vars).astype(np.float32))
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, var=var)
    input_path = os.path.join(temp_dir, "input.h5ad")
    adata.write_h5ad(input_path)

    # Create gene lists
    genelists = {"test_list": ["gene_0", "gene_1"]}
    genelists_path = os.path.join(temp_dir, "genelists.json")
    with open(genelists_path, "w") as f:
        json.dump(genelists, f)

    output_path = os.path.join(temp_dir, "output.h5ad")

    result = subprocess.run(
        [
            "adata_util",
            "score_genes",
            input_path,
            genelists_path,
            output_path,
            "--ctrl_size",
            "10",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "ctrl_size" in result.stdout

    adata_out = ad.read_h5ad(output_path)
    assert "test_list_score" in adata_out.obs.columns
