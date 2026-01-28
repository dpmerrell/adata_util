"""Scanpy command: run scanpy operations on AnnData objects."""

import scanpy as sc
import anndata as ad


def register(subparsers):
    """Register the scanpy subcommand."""
    parser = subparsers.add_parser(
        "scanpy",
        help="Run a scanpy operation on an AnnData in an h5ad file",
    )
    parser.add_argument(
        "operation",
        help="Scanpy operation to run (e.g., 'pp.normalize_total', 'tl.pca')",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        help="Path to the output h5ad file. If not specified, overwrites input.",
    )
    parser.set_defaults(func=run)


def _parse_value(value_str):
    """Parse a string value into appropriate Python type."""
    if value_str.lower() == "true":
        return True
    if value_str.lower() == "false":
        return False
    if value_str.lower() == "none":
        return None

    try:
        return int(value_str)
    except ValueError:
        pass

    try:
        return float(value_str)
    except ValueError:
        pass

    return value_str


def _parse_extra_args(extra_args):
    """Parse a list of extra arguments into kwargs dict."""
    kwargs = {}
    i = 0
    while i < len(extra_args):
        arg = extra_args[i]
        if arg.startswith("--"):
            key = arg[2:].replace("-", "_")
            if i + 1 < len(extra_args) and not extra_args[i + 1].startswith("--"):
                kwargs[key] = _parse_value(extra_args[i + 1])
                i += 2
            else:
                kwargs[key] = True
                i += 1
        else:
            i += 1
    return kwargs


def _get_scanpy_function(operation):
    """Get a scanpy function from its dotted name."""
    parts = operation.split(".")
    obj = sc
    for part in parts:
        obj = getattr(obj, part)
    return obj


def run(args, extra_args=None):
    """Execute the scanpy command."""
    kwargs = _parse_extra_args(extra_args or [])

    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    print(f"Running scanpy.{args.operation}...")
    try:
        func = _get_scanpy_function(args.operation)
    except AttributeError:
        print(f"Error: Unknown scanpy operation '{args.operation}'")
        return

    if kwargs:
        print(f"  with kwargs: {kwargs}")
    result = func(adata, **kwargs)

    if result is not None and isinstance(result, ad.AnnData):
        adata = result

    output = args.output if args.output else args.input

    print(f"Writing to {output}...")
    adata.write_h5ad(output)
    print("Done.")
