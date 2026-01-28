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
    # Handle boolean
    if value_str.lower() == "true":
        return True
    if value_str.lower() == "false":
        return False
    if value_str.lower() == "none":
        return None

    # Try integer
    try:
        return int(value_str)
    except ValueError:
        pass

    # Try float
    try:
        return float(value_str)
    except ValueError:
        pass

    # Return as string
    return value_str


def _get_scanpy_function(operation):
    """Get a scanpy function from its dotted name."""
    parts = operation.split(".")
    obj = sc
    for part in parts:
        obj = getattr(obj, part)
    return obj


def run(args, extra_args=None):
    """Execute the scanpy command."""
    # Parse extra arguments as kwargs
    # We need to use parse_known_args to capture extra --key value pairs
    import sys

    # Get the full argv and find our operation
    argv = sys.argv
    try:
        op_idx = argv.index(args.operation)
    except ValueError:
        op_idx = 0

    # Everything after input/output files are kwargs
    # Find where the known args end
    kwargs = {}
    i = 0
    remaining = argv[op_idx + 1:]  # After operation name

    # Skip input and output files
    file_count = 0
    while i < len(remaining) and not remaining[i].startswith("--"):
        file_count += 1
        i += 1

    # Parse remaining as kwargs
    while i < len(remaining):
        arg = remaining[i]
        if arg.startswith("--"):
            key = arg[2:].replace("-", "_")
            if i + 1 < len(remaining) and not remaining[i + 1].startswith("--"):
                value = _parse_value(remaining[i + 1])
                i += 2
            else:
                # Flag without value, treat as True
                value = True
                i += 1
            kwargs[key] = value
        else:
            i += 1

    # Read input
    print(f"Reading {args.input}...")
    adata = ad.read_h5ad(args.input)

    # Get the scanpy function
    print(f"Running scanpy.{args.operation}...")
    try:
        func = _get_scanpy_function(args.operation)
    except AttributeError:
        print(f"Error: Unknown scanpy operation '{args.operation}'")
        return

    # Call the function
    if kwargs:
        print(f"  with kwargs: {kwargs}")
    result = func(adata, **kwargs)

    # Some functions return a modified adata, others modify in place
    if result is not None and isinstance(result, ad.AnnData):
        adata = result

    # Determine output path
    output = args.output if args.output else args.input

    # Write output
    print(f"Writing to {output}...")
    adata.write_h5ad(output)
    print("Done.")
