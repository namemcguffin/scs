import sys
import os
from pathlib import Path
from anndata import read_h5ad, AnnData
from pandas import DataFrame

def from_anndata(adata: AnnData, out_path: Path, spatial_slot: str = "spatial", layer: str | None = None):
    spatial_data = adata.obsm[spatial_slot]
    mtx = adata.layers[layer] if layer is not None else adata.X

    os.makedirs(out_path, exist_ok=True)
    os.makedirs(out_path / "feat", exist_ok=True)

    DataFrame(
        {"cell": adata.obs_names, "x": spatial_data[:, 0], "y": spatial_data[:, 1]}
    ).to_csv(out_path / "cells.tsv", sep="\t", index=False)

    for i, feat in enumerate(adata.var_names):
        nz = mtx[:, i].nonzero()[0]
        DataFrame(
            {
                "cell": (adata.obs_names[j] for j in nz),
                "expr": (mtx[j, i] for j in nz),
            }
        ).to_csv(os.path.join(out_path, "feat", f"{feat}.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    from_anndata(
        read_h5ad(sys.argv[1]),
        Path(sys.argv[2]),
        *(sys.argv[3:5])
    )
