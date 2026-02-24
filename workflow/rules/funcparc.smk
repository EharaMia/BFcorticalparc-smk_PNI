from pathlib import Path

PNI_ROOT = "/home/m-ehara/project/PNI/derivatives"

# BF マスクをコピー
rule copy_bf_seed:
    input:
        seed="/home/m-ehara/project/PNI/resource/seed_1p6mm.nii.gz"
    output:
        bf_mask="results/funcparc/atlas/tpl-MNI152NLin6Asym_label-fullBF_mask.nii.gz"
    shell:
        """
        mkdir -p results/funcparc/atlas
        cp {input.seed} {output.bf_mask}
        """

# BF voxel 時系列（599 × T）
rule prepare_subcortical:
    input:
        fmri=(
            PNI_ROOT
            + "/sub-{subject}/ses-{session}/func/desc-me_task-rest_bold/volumetric/"
            + "sub-{subject}_ses-{session}_space-func_desc-me_preproc.nii.gz"
        ),
        bf_mask="results/funcparc/atlas/tpl-MNI152NLin6Asym_label-fullBF_mask.nii.gz"
    output:
        bf_ts="results/funcparc/sub-{subject}/ses-{session}/BF_voxel_timeseries.npy"
    script:
        "../scripts/prepare_subcortical.py"

# cortex parcel 時系列（400 × T）
rule make_cortex_parcel_ts:
    input:
        ts_gii=(
            PNI_ROOT
            + "/sub-{subject}/ses-{session}/func/desc-me_task-rest_bold/surf/"
            + "sub-{subject}_ses-{session}_surf-fsLR-32k_desc-timeseries_clean.shape.gii"
        ),
        lab_L="/home/m-ehara/project/PNI/resource/schaefer_32k/Schaefer2018_7Networks_400.32k.L.label.gii",
        lab_R="/home/m-ehara/project/PNI/resource/schaefer_32k/Schaefer2018_7Networks_400.32k.R.label.gii",
    output:
        cortex_ts="results/funcparc/sub-{subject}/ses-{session}/cortex_parcel_timeseries.npy"
    script:
        "../scripts/make_cortex_parcel_ts.py"


# BF × cortex 相関（599 × 400）
rule compute_correlation:
    input:
        bf_ts="results/funcparc/sub-{subject}/ses-{session}/BF_voxel_timeseries.npy",
        cortex_ts="results/funcparc/sub-{subject}/ses-02/cortex_parcel_timeseries.npy"
    output:
        corr="results/funcparc/sub-{subject}/ses-{session}/"
             "sub-{subject}_ses-{session}_seed-BF_schaefer400_corr_voxelwise.npz"
    script:
        "../scripts/compute_correlation.py"
  
# 10 人分まとめる
rule combine_correlation:
    """Concatenate correlation matrices across subject"""
    output:
        combined_corr="results/funcparc/group/clustering/"
                      "tpl-MNI152NLin6Asym_label-fullBF_desc-correlationMatrix_1p6mm.npz"
    params:
        pattern="results/funcparc/sub-*/ses-02/sub-*_ses-02_seed-BF_schaefer400_corr_voxelwise.npz"
    container:
        config["singularity"]["pythondeps"]
    group:
        "funcparc_group2"
    script:
        "../scripts/combine_correlation.py"

