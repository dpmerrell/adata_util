"""Command line interface for adata_util."""

import argparse
import sys

from adata_util.commands import view, concat, join, scanpy_cmd


def main():
    """Main entry point for the adata_util CLI."""
    parser = argparse.ArgumentParser(
        prog="adata_util",
        description="A CLI for operating on AnnData objects (h5ad files).",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Register subcommands
    view.register(subparsers)
    concat.register(subparsers)
    join.register(subparsers)
    scanpy_cmd.register(subparsers)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Execute the command
    args.func(args)


if __name__ == "__main__":
    main()
