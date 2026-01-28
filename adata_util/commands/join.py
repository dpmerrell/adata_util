"""Join command: join tabular data to AnnData obs or var."""

import anndata as ad
import pandas as pd


def register(subparsers):
    """Register the join subcommand."""
    parser = subparsers.add_parser(
        "join",
        help="Join obs- or var-level data from a tabular file to an AnnData",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.add_argument(
        "table",
        help="Path to the tabular file (csv, tsv, or tab-delimited)",
    )
    parser.add_argument(
        "output",
        help="Path to the output h5ad file",
    )
    parser.add_argument(
        "--on",
        default="obs",
        choices=["obs", "var"],
        help="Join to obs (cell metadata) or var (gene metadata). Default: obs",
    )
    parser.add_argument(
        "--left-on",
        default=None,
        help="Column in AnnData obs/var to join on. Default: index",
    )
    parser.add_argument(
        "--right-on",
        default=None,
        help="Column in table to join on. Default: first column",
    )
    parser.add_argument(
        "--how",
        default="left",
        choices=["left", "inner", "outer"],
        help="Type of join. Default: left",
    )
    parser.add_argument(
        "--sep",
        default=None,
        help="Delimiter for table file. Default: auto-detect (comma or tab)",
    )
    parser.set_defaults(func=run)


def _read_table(path, sep=None):
    """Read a tabular file, auto-detecting delimiter if not specified."""
    if sep is not None:
        return pd.read_csv(path, sep=sep)

    # Auto-detect based on extension
    if path.endswith(".csv"):
        return pd.read_csv(path)
    elif path.endswith(".tsv") or path.endswith(".tab"):
        return pd.read_csv(path, sep="\t")
    else:
        # Try to detect from first line
        with open(path, "r") as f:
            first_line = f.readline()
        if "\t" in first_line:
            return pd.read_csv(path, sep="\t")
        else:
            return pd.read_csv(path)


def run(args):
    """Execute the join command."""
    # Read input files
    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    print(f"Reading {args.table}...")
    table = _read_table(args.table, args.sep)

    # Determine join keys
    right_on = args.right_on
    if right_on is None:
        right_on = table.columns[0]

    # Set the right key as index for the table
    table = table.set_index(right_on)

    # Get the target dataframe (obs or var)
    if args.on == "obs":
        target = adata.obs
    else:
        target = adata.var

    # Determine left key
    left_on = args.left_on

    if left_on is not None:
        # Use the specified column as join key
        original_index = target.index.copy()
        target = target.reset_index()
        target = target.set_index(left_on)

    # Perform the join
    print(f"Joining on {args.on}...")

    # Find columns to add (avoid duplicates)
    new_cols = [c for c in table.columns if c not in target.columns]
    if len(new_cols) == 0:
        print("Warning: No new columns to add (all columns already exist)")
    else:
        print(f"Adding columns: {', '.join(new_cols)}")

    # Join the tables
    result = target.join(table[new_cols], how=args.how)

    # Restore original index if we changed it
    if left_on is not None:
        result = result.reset_index()
        result = result.set_index(original_index.name if original_index.name else "index")
        result.index = original_index

    # Update the AnnData
    if args.on == "obs":
        adata.obs = result
    else:
        adata.var = result

    # Write output
    print(f"Writing to {args.output}...")
    adata.write_h5ad(args.output)
    print("Done.")
