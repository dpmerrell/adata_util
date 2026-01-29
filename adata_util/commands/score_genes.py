"""Score genes command: score gene lists using scanpy.tl.score_genes."""

import json

import anndata as ad
import scanpy as sc

from adata_util.commands.utils import parse_extra_args


def register(subparsers):
    """Register the score_genes subcommand."""
    parser = subparsers.add_parser(
        "score_genes",
        help="Score gene lists using scanpy.tl.score_genes",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.add_argument(
        "genelists",
        help="Path to a YAML or JSON file containing gene lists",
    )
    parser.add_argument(
        "output",
        help="Path to the output h5ad file",
    )
    parser.set_defaults(func=run)


def _load_genelists(path):
    """Load gene lists from a YAML or JSON file."""
    if path.endswith(".yaml") or path.endswith(".yml"):
        import yaml
        with open(path, "r") as f:
            return yaml.safe_load(f)
    else:
        with open(path, "r") as f:
            return json.load(f)


def run(args, extra_args=None):
    """Execute the score_genes command."""
    kwargs = parse_extra_args(extra_args or [])

    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    print(f"Loading gene lists from {args.genelists}...")
    genelists = _load_genelists(args.genelists)

    if not isinstance(genelists, dict):
        print("Error: genelists file must contain a dictionary")
        return

    print(f"Scoring {len(genelists)} gene lists...")
    if kwargs:
        print(f"  with kwargs: {kwargs}")

    for name, genes in genelists.items():
        score_name = f"{name}_score"
        print(f"  {name}: {len(genes)} genes -> {score_name}")
        sc.tl.score_genes(adata, genes, score_name=score_name, **kwargs)

    print(f"Writing to {args.output}...")
    adata.write_h5ad(args.output)
    print("Done.")
