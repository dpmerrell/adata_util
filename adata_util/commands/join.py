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
        choices=["left", "inner"],
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


def run(args, _extra_args=None):
    """Execute the join command."""
    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    print(f"Reading {args.table}...")
    table = _read_table(args.table, args.sep)

    # Determine join keys
    right_on = args.right_on if args.right_on is not None else table.columns[0]
    table = table.set_index(right_on)

    # Get target dataframe and set up join key
    target = adata.obs.copy() if args.on == "obs" else adata.var.copy()
    left_on = args.left_on

    if left_on is not None:
        target = target.reset_index().set_index(left_on)

    # Perform the join
    print(f"Joining on {args.on} ({args.how})...")
    new_cols = [c for c in table.columns if c not in target.columns]
    if new_cols:
        print(f"Adding columns: {', '.join(new_cols)}")
    else:
        print("Warning: No new columns to add (all columns already exist)")

    result = target.join(table[new_cols] if new_cols else table[[]], how=args.how)

    # Subset AnnData to keys in the result
    keys_to_keep = result.index
    if args.on == "obs":
        if left_on is not None:
            mask = adata.obs[left_on].isin(keys_to_keep)
        else:
            mask = adata.obs_names.isin(keys_to_keep)
        adata = adata[mask, :].copy()
    else:
        if left_on is not None:
            mask = adata.var[left_on].isin(keys_to_keep)
        else:
            mask = adata.var_names.isin(keys_to_keep)
        adata = adata[:, mask].copy()

    # Restore original index structure and assign to AnnData
    if left_on is not None:
        result = result.reset_index()
        if args.on == "obs":
            result.index = adata.obs_names
        else:
            result.index = adata.var_names

    if args.on == "obs":
        adata.obs = result
    else:
        adata.var = result

    print(f"Writing to {args.output}...")
    adata.write_h5ad(args.output)
    print("Done.")
