"""Plot embedding command: plot embeddings from AnnData using scanpy."""

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from adata_util.commands.utils import parse_extra_args


def register(subparsers):
    """Register the plot_embedding subcommand."""
    parser = subparsers.add_parser(
        "plot_embedding",
        help="Plot an embedding from an h5ad file using scanpy.pl.embedding",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.add_argument(
        "basis",
        help="Basis for the embedding (e.g., 'pca', 'umap', 'tsne')",
    )
    parser.add_argument(
        "output",
        help="Path to the output image file (e.g., .png, .pdf)",
    )
    parser.set_defaults(func=run)


def run(args, extra_args=None):
    """Execute the plot_embedding command."""
    kwargs = parse_extra_args(extra_args or [])

    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    print(f"Plotting embedding (basis={args.basis})...")
    if kwargs:
        print(f"  with kwargs: {kwargs}")

    sc.pl.embedding(adata, basis=args.basis, show=False, **kwargs)

    print(f"Saving to {args.output}...")
    plt.savefig(args.output, bbox_inches="tight")
    plt.close()
    print("Done.")
