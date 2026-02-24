import numpy as np

inp = snakemake.input.group_npz
out = snakemake.output.group_npz_400

DROP_IDXS = [0, 201]   # label-1000, label-2000

d = np.load(inp)
conn = d["conn_group"]   # (10,599,402)
mask = d["mask"]
aff  = d["affine"]

keep = np.ones(conn.shape[2], dtype=bool)
keep[DROP_IDXS] = False

conn400 = conn[:, :, keep]

with open(snakemake.log[0], "w") as f:
    f.write(f"input shape: {conn.shape}\n")
    f.write(f"drop idxs: {DROP_IDXS}\n")
    f.write(f"output shape: {conn400.shape}\n")

np.savez(out, conn_group=conn400, mask=mask, affine=aff)
