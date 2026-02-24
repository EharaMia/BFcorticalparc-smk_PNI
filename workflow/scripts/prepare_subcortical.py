import nibabel as nib
import numpy as np

# BF とrsfMRI を用意
fmri = nib.load(snakemake.input["fmri"])   # rsfMRI の4次元画像（x, y, z, time）
bf   = nib.load(snakemake.input["bf_mask"])   # BFマスク

# 配列
fmri_data = fmri.get_fdata()
bf_data   = bf.get_fdata() > 0   # BF のみ取り出す

# 空間定義
fmri_aff = fmri.affine
bf_aff   = bf.affine

# BF ボクセルの場所を取り出す
bf_vox = np.argwhere(bf_data)

# BF ボクセルをmm に変換
bf_xyz = nib.affines.apply_affine(bf_aff, bf_vox)

# rsfMRI 画像のボクセル座標との位置対応を合わせる
fmri_ijk = nib.affines.apply_affine(
    np.linalg.inv(fmri_aff), bf_xyz
)
fmri_ijk = np.round(fmri_ijk).astype(int)

# 時系列の箱を作る（行：BFボクセル(599), 列：時間点(T), 初期値：NaN）
T = fmri_data.shape[3]
ts = np.full((bf_vox.shape[0], T), np.nan)

# BOLD 時系列を拾う
for i, (x, y, z) in enumerate(fmri_ijk):   # BF599 ボクセルを1 つずつ見る
    # refMRI の範囲内ならBOLD 時系列をコピー、範囲外な NaN のまま
    if (
        0 <= x < fmri_data.shape[0] and
        0 <= y < fmri_data.shape[1] and
        0 <= z < fmri_data.shape[2]
    ):
        ts[i] = fmri_data[x, y, z, :]

# 保存
np.save(snakemake.output["bf_ts"], ts)
