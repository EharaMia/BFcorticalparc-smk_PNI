import numpy as np

inp = snakemake.input.group_npz
out = snakemake.output.group_npz_400

d = np.load(inp)
conn = d["conn_group"]   # (nsub, 599, 402)
mask = d["mask"]
aff  = d["affine"]

# 全被験者×全voxelで完全ゼロの列を探す
is_all_zero = (conn == 0).all(axis=(0, 1))   # shape=(402,)
zero_idx = np.where(is_all_zero)[0]

with open(snakemake.log[0], "w") as f:
    f.write(f"input shape: {conn.shape}\n")
    f.write(f"all-zero target indices: {zero_idx.tolist()}\n")
    f.write(f"n all-zero targets: {len(zero_idx)}\n")

conn400 = conn[:, :, ~is_all_zero]

with open(snakemake.log[0], "a") as f:
    f.write(f"output shape: {conn400.shape}\n")

np.savez(out, conn_group=conn400, mask=mask, affine=aff)
