"""View command: prints an informative summary of an h5ad file."""

import anndata as ad


def register(subparsers):
    """Register the view subcommand."""
    parser = subparsers.add_parser(
        "view",
        help="Print an informative summary of the contents of an h5ad file",
    )
    parser.add_argument(
        "input",
        help="Path to the input h5ad file",
    )
    parser.set_defaults(func=run)


def run(args):
    """Execute the view command."""
    adata = ad.read_h5ad(args.input)

    print(f"AnnData object: {args.input}")
    print(f"  Shape: {adata.n_obs} observations x {adata.n_vars} variables")
    print()

    # obs (cell metadata)
    print(f"obs (observation/cell metadata): {len(adata.obs.columns)} columns")
    if len(adata.obs.columns) > 0:
        for col in adata.obs.columns:
            dtype = adata.obs[col].dtype
            print(f"    {col}: {dtype}")
    print()

    # var (gene metadata)
    print(f"var (variable/gene metadata): {len(adata.var.columns)} columns")
    if len(adata.var.columns) > 0:
        for col in adata.var.columns:
            dtype = adata.var[col].dtype
            print(f"    {col}: {dtype}")
    print()

    # obsm (observation matrices)
    print(f"obsm (observation matrices): {len(adata.obsm.keys())} entries")
    if len(adata.obsm.keys()) > 0:
        for key in adata.obsm.keys():
            shape = adata.obsm[key].shape
            print(f"    {key}: {shape}")
    print()

    # varm (variable matrices)
    print(f"varm (variable matrices): {len(adata.varm.keys())} entries")
    if len(adata.varm.keys()) > 0:
        for key in adata.varm.keys():
            shape = adata.varm[key].shape
            print(f"    {key}: {shape}")
    print()

    # obsp (observation pairwise)
    print(f"obsp (observation pairwise): {len(adata.obsp.keys())} entries")
    if len(adata.obsp.keys()) > 0:
        for key in adata.obsp.keys():
            shape = adata.obsp[key].shape
            print(f"    {key}: {shape}")
    print()

    # varp (variable pairwise)
    print(f"varp (variable pairwise): {len(adata.varp.keys())} entries")
    if len(adata.varp.keys()) > 0:
        for key in adata.varp.keys():
            shape = adata.varp[key].shape
            print(f"    {key}: {shape}")
    print()

    # uns (unstructured)
    print(f"uns (unstructured): {len(adata.uns.keys())} entries")
    if len(adata.uns.keys()) > 0:
        for key in adata.uns.keys():
            val = adata.uns[key]
            val_type = type(val).__name__
            print(f"    {key}: {val_type}")
    print()

    # layers
    print(f"layers: {len(adata.layers.keys())} entries")
    if len(adata.layers.keys()) > 0:
        for key in adata.layers.keys():
            shape = adata.layers[key].shape
            print(f"    {key}: {shape}")
