"""Plot command: run scanpy plotting functions."""

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from adata_util.commands.utils import parse_extra_args


def register(subparsers):
    """Register the plot subcommand."""
    parser = subparsers.add_parser(
        "plot",
        help="Run a scanpy.pl plotting function on an h5ad file",
    )
    parser.add_argument(
        "function",
        help="Plotting function to run (e.g., 'embedding', 'violin', 'heatmap')",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.add_argument(
        "output",
        help="Path to the output image file (e.g., .png, .pdf)",
    )
    parser.set_defaults(func=run)


def _get_plot_function(function_name):
    """Get a scanpy.pl function from its name."""
    parts = function_name.split(".")
    obj = sc.pl
    for part in parts:
        obj = getattr(obj, part)
    return obj


def run(args, extra_args=None):
    """Execute the plot command."""
    kwargs = parse_extra_args(extra_args or [])

    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    print(f"Running scanpy.pl.{args.function}...")
    try:
        func = _get_plot_function(args.function)
    except AttributeError:
        print(f"Error: Unknown plotting function 'scanpy.pl.{args.function}'")
        return

    if kwargs:
        print(f"  with kwargs: {kwargs}")

    func(adata, show=False, **kwargs)

    print(f"Saving to {args.output}...")
    plt.savefig(args.output, bbox_inches="tight")
    plt.close()
    print("Done.")
