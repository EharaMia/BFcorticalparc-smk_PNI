import numpy as np

bf_ts = np.load(snakemake.input[0])        # (599, T)
ctx_ts = np.load(snakemake.input[1])       # (400, T or F)

# 次元を揃える（最小長）
T = min(bf_ts.shape[1], ctx_ts.shape[1])
bf_ts  = bf_ts[:, :T]
ctx_ts = ctx_ts[:, :T]

# 平均除去
bf_ts  -= bf_ts.mean(axis=1, keepdims=True)
ctx_ts -= ctx_ts.mean(axis=1, keepdims=True)

# ノルム
bf_norm  = np.linalg.norm(bf_ts, axis=1, keepdims=True)
ctx_norm = np.linalg.norm(ctx_ts, axis=1, keepdims=True)

# 相関
corr = (bf_ts @ ctx_ts.T) / (bf_norm * ctx_norm.T)

# 数値安全
corr[~np.isfinite(corr)] = 0.0

# Fisher z
corr_z = np.arctanh(np.clip(corr, -0.999999, 0.999999))

np.savez(
    snakemake.output[0],
    corr=corr_z,
    n_bf_voxels=bf_ts.shape[0],
    n_parcels=ctx_ts.shape[0],
    n_features=T,
)
