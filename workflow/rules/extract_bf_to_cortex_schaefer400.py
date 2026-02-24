#!/usr/bin/env python3
import numpy as np


def main(npz_path, output_txt):
    # ===== load =====
    d = np.load(npz_path)
    corr = d["corr"]

    # ===== safety checks =====
    if corr.ndim != 3:
        raise ValueError("corr must be (n_subjects, N, N)")
    if corr.shape[0] != 1:
        raise ValueError("currently expecting one subject")

    N = corr.shape[1]
    if N < 401:
        raise ValueError("Matrix too small to contain BF + Schaefer-400")

    # ===== indices (HCP-compatible) =====
    bf_idx = 0
    cortex_start = N - 400
    cortex_end = N

    # ===== extract =====
    bf_to_cortex = corr[0, bf_idx, cortex_start:cortex_end]

    if bf_to_cortex.shape != (400,):
        raise ValueError("Expected shape (400,)")

    # ===== save =====
    np.savetxt(output_txt, bf_to_cortex[None, :], fmt="%.6f")

    print("✅ BF → cortex extracted:", bf_to_cortex.shape)
    print("range:", bf_to_cortex.min(), "to", bf_to_cortex.max())


if __name__ == "__main__":
    main(
        npz_path=snakemake.input.corr,
        output_txt=snakemake.output.vec,
    )
