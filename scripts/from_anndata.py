import sys
import os
from anndata import read_h5ad
from pandas import DataFrame

if __name__ == "__main__":
    adata = read_h5ad(sys.argv[1])
    out_path = sys.argv[2]
    spatial_data = adata.obsm[sys.argv[3] if len(sys.argv) >= 4 else "spatial"]
    mtx = adata.layers[sys.argv[4]] if len(sys.argv) >= 5 else adata.X

    os.makedirs(out_path, exist_ok=True)
    os.makedirs(os.path.join(out_path, "feat"), exist_ok=True)

    DataFrame(
        {"cell": adata.obs_names, "x": spatial_data[:, 0], "y": spatial_data[:, 1]}
    ).to_csv(os.path.join(out_path, "cells.tsv"), sep="\t", index=False)

    for i, feat in enumerate(adata.var_names):
        nz = mtx[:, i].nonzero()[0]
        DataFrame(
            {
                "cell": (adata.obs_names[j] for j in nz),
                "expr": (mtx[j, i] for j in nz),
            }
        ).to_csv(os.path.join(out_path, "feat", f"{feat}.tsv"), sep="\t", index=False)
