import numpy as np
import nibabel as nib

# timeseries_clean: (T, 64984)
g = nib.load(snakemake.input["ts_gii"])
ts = g.darrays[0].data

if ts.ndim != 2:
    raise RuntimeError(f"Unexpected ts dims: {ts.shape}")

T, V = ts.shape
print("timeseries shape:", ts.shape)

# label: (32492,)
labL = nib.load(snakemake.input["lab_L"]).darrays[0].data.astype(int)
labR = nib.load(snakemake.input["lab_R"]).darrays[0].data.astype(int)

if labL.ndim != 1 or labR.ndim != 1:
    raise RuntimeError(f"Unexpected label dims: L={labL.shape} R={labR.shape}")

Vh = labL.shape[0]
if labR.shape[0] != Vh:
    raise RuntimeError("L/R label vertex counts differ")

if V != 2 * Vh:
    raise RuntimeError(f"Vertex mismatch: ts has {V}, labels imply {2*Vh}")

# 64984頂点のラベル配列を作る
# Left: 1..200
# Right: 1..200 を 201..400 にずらす
labels = np.zeros((V,), dtype=int)
labels[:Vh] = labL
labels[Vh:] = labR

# (V, T) にして平均しやすくする
tsV = ts.T  # (64984, 210)

out = np.full((400, T), np.nan, dtype=np.float32)

for k in range(1, 401):
    m = (labels == k)
    if not np.any(m):
        continue
    out[k-1] = tsV[m].mean(axis=0)

missing = np.where(~np.isfinite(out).all(axis=1))[0]
print("missing parcels:", len(missing), "example:", missing[:20].tolist())

if len(missing) > 0:
    raise RuntimeError(f"Some parcels missing: count={len(missing)} example={missing[:20].tolist()}")

np.save(snakemake.output["cortex_ts"], out)
print("Saved cortex_ts:", out.shape)
