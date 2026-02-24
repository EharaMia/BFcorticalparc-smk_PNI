import nibabel as nib
import numpy as np

# 1) seed mask 読み込み
mask_nib = nib.load(snakemake.input.mask)
mask_vol = mask_nib.get_fdata()
mask_indices = mask_vol > 0

nvoxels = int(mask_indices.sum())

# 2) targets.txt から connmap ファイル一覧を読む
with open(snakemake.params.targets_txt, "r") as f:
    target_paths = [line.strip() for line in f if line.strip()]

# probtrackx出力のprefixに合わせて置換
# targets.txt は ".../targets/sub-PNCxxx_label-1000_mask.nii.gz"
# Linear出力は ".../probtrackx_output_seed-.../seeds_to_sub-PNCxxx_label-1000_mask.nii.gz"
conn_files = []
for t in target_paths:
    base = t.split("/")[-1]  # sub-PNCxxx_label-1000_mask.nii.gz
    conn_files.append(f"{snakemake.params.probtrack_dir}/seeds_to_{base}")

ntargets = len(conn_files)

# 3) conn matrix 作成
conn = np.zeros((nvoxels, ntargets), dtype=np.float32)

for i, conn_file in enumerate(conn_files):
    vol = nib.load(conn_file).get_fdata()
    conn[:, i] = vol[mask_indices]

np.savez(
    snakemake.output.connmap_npz,
    conn=conn,
    mask=mask_vol,
    affine=mask_nib.affine,
)
