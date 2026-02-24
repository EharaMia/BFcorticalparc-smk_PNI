#!/usr/bin/env python3
import glob
import numpy as np

pattern = snakemake.params.pattern
out = snakemake.output.combined_corr

files = sorted(glob.glob(pattern))
if len(files) == 0:
    raise FileNotFoundError(f"No files matched pattern: {pattern}")

first = np.load(files[0])
if "corr" not in first.files:
    raise KeyError(f"{files[0]} does not contain 'corr'. keys={first.files}")

shape = first["corr"].shape
combined = np.zeros((len(files), shape[0], shape[1]), dtype=np.float32)

for i, f in enumerate(files):
    d = np.load(f)
    if d["corr"].shape != shape:
        raise ValueError(f"Shape mismatch: {f} has {d['corr'].shape}, expected {shape}")
    combined[i] = d["corr"]

# メタ情報は「存在するものだけ」保存する
meta = {}
for k in ["n_bf_voxels", "n_parcels", "n_features", "n_timepoints"]:
    if k in first.files:
        meta[k] = first[k]

np.savez(out, corr_group=combined, **meta)

