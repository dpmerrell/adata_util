"""Split command: extract obs or var table from AnnData to a tabular file."""

import anndata as ad

from adata_util.commands.utils import parse_extra_args


def register(subparsers):
    """Register the split subcommand."""
    parser = subparsers.add_parser(
        "split",
        help="Write obs or var table from an h5ad file to a tabular file (CSV/TSV)",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.add_argument(
        "table",
        choices=["obs", "var"],
        help="Which table to extract (obs or var)",
    )
    parser.add_argument(
        "output",
        help="Path to the output tabular file (.csv or .tsv)",
    )
    parser.add_argument(
        "--kept-cols",
        nargs="+",
        default=None,
        help="Columns to keep (default: all columns). Index is always included.",
    )
    parser.set_defaults(func=run)


def run(args, _extra_args=None):
    """Execute the split command."""
    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    # Get the target table
    if args.table == "obs":
        df = adata.obs.copy()
    else:
        df = adata.var.copy()

    # Filter columns if specified
    if args.kept_cols is not None:
        missing = [c for c in args.kept_cols if c not in df.columns]
        if missing:
            print(f"Error: columns not found: {', '.join(missing)}")
            return
        df = df[args.kept_cols]
        print(f"Keeping columns: {', '.join(args.kept_cols)}")
    else:
        print(f"Writing all {len(df.columns)} columns")

    # Determine output format based on extension
    if args.output.endswith(".tsv") or args.output.endswith(".tab"):
        sep = "\t"
    else:
        sep = ","

    # Write output (index is always included)
    print(f"Writing to {args.output}...")
    df.to_csv(args.output, sep=sep, index=True)
    print("Done.")
