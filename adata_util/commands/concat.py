"""Concat command: concatenates AnnData objects from multiple h5ad files."""

import anndata as ad


def register(subparsers):
    """Register the concat subcommand."""
    parser = subparsers.add_parser(
        "concat",
        help="Concatenate AnnData objects from one or more h5ad files",
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="Input h5ad files followed by output h5ad file",
    )
    parser.add_argument(
        "--axis",
        type=int,
        default=0,
        choices=[0, 1],
        help="Axis along which to concatenate (0=obs, 1=var). Default: 0",
    )
    parser.add_argument(
        "--join",
        default="inner",
        choices=["inner", "outer"],
        help="How to handle the other axis. Default: inner",
    )
    parser.add_argument(
        "--label",
        default=None,
        help="Column name for batch labels in obs/var. Default: None",
    )
    parser.set_defaults(func=run)


def run(args):
    """Execute the concat command."""
    if len(args.files) < 2:
        print("Error: concat requires at least one input file and one output file")
        return

    # Last argument is the output file
    input_files = args.files[:-1]
    output_file = args.files[-1]

    # Read all input files
    adatas = []
    for f in input_files:
        print(f"Reading {f}...")
        adatas.append(ad.read_h5ad(f))

    # Concatenate
    print(f"Concatenating {len(adatas)} AnnData objects...")

    concat_kwargs = {
        "axis": args.axis,
        "join": args.join,
    }
    if args.label is not None:
        concat_kwargs["label"] = args.label

    result = ad.concat(adatas, **concat_kwargs)

    # Write output
    print(f"Writing to {output_file}...")
    result.write_h5ad(output_file)
    print("Done.")
