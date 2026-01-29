"""Command line interface for adata_util."""

import argparse
import sys

import argcomplete

from adata_util.commands import view, concat, join, scanpy_cmd, plot_embedding, split, score_genes


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
    plot_embedding.register(subparsers)
    split.register(subparsers)
    score_genes.register(subparsers)

    argcomplete.autocomplete(parser)
    args, remaining = parser.parse_known_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Execute the command, passing remaining args for commands that accept them
    args.func(args, remaining)


if __name__ == "__main__":
    main()

