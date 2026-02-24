rule drop_1000_2000_targets_group_pni:
    input:
        group_npz = "results/diffparc/tpl-{template}/tpl-{template}_label-fullBF_desc-Linear_schaefer400_from-group_conn_voxelwise.npz"
    output:
        group_npz_400 = "results/diffparc/tpl-{template}/tpl-{template}_label-fullBF_desc-Linear_schaefer400_from-group_conn_voxelwise_400targets.npz"
    log:
        "logs/drop_1000_2000_targets_group/fullBF_{template}.log"
    conda:
        "../envs/sklearn.yml"
    script:
        "../scripts/drop_1000_2000_targets_group_pni.py"


rule gather_connmap_group_pni:
    input:
        connmap_npz = expand(
            "results/diffparc/sub-{subject}/sub-{subject}_seed-fullBF_space-{template}_desc-Linear_schaefer400_conn_voxelwise.npz",
            subject=subjects,
            template="{template}"
        )
    output:
        connmap_group_npz = "results/diffparc/tpl-{template}/tpl-{template}_label-fullBF_desc-Linear_schaefer400_from-group_conn_voxelwise.npz"
    log:
        "logs/gather_connmap_group/fullBF_{template}_Linear_schaefer400.log"
    conda:
        "../envs/sklearn.yml"
    group:
        "group1"
    script:
        "../scripts/gather_connmap_group_pni.py"


rule save_connmap_template_npz_linear:
    input:
        mask = bids(root='results/diffparc', template='{template}', label='{seed}', suffix='mask.nii.gz'),
        done = "results/diffparc/sub-{subject}/probtrackx_output_seed-{seed}_space-{template}_desc-Linear/DONE.txt"
    params:
        probtrack_dir = "results/diffparc/sub-{subject}/probtrackx_output_seed-{seed}_space-{template}_desc-Linear",
        targets_txt = "results/diffparc/sub-{subject}/targets.txt"
    output:
        connmap_npz = "results/diffparc/sub-{subject}/sub-{subject}_seed-{seed}_space-{template}_desc-Linear_schaefer400_conn_voxelwise.npz"
    log:
        "logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}_Linear.log"
    group: "participant1"
    conda: "../envs/sklearn.yml"
    script: "../scripts/save_connmap_template_npz_pni.py"





from snakemake.io import glob_wildcards

rule transform_conn_to_template_linear:
    input:
        probtrack_dir = "results/diffparc/sub-{subject}/probtrackx_output",
        warp = config['ants_warp_nii'],
        ref = bids(root='results/diffparc', template='{template}', label='{seed}', suffix='mask.nii.gz')
    params:
        in_files = lambda wildcards, input: expand(
            input.probtrack_dir + "/seeds_to_sub-" + wildcards.subject + "_label-{lab}_mask.nii.gz",
            lab=sorted(
                glob_wildcards(
                    input.probtrack_dir + "/seeds_to_sub-" + wildcards.subject + "_label-{lab}_mask.nii.gz"
                ).lab
            )
        ),
        out_files = lambda wildcards, output: expand(
            "results/diffparc/sub-{subject}/probtrackx_output_seed-{seed}_space-{template}_desc-Linear/"
            "seeds_to_sub-{subject}_label-{lab}_mask.nii.gz",
            subject=wildcards.subject,
            seed=wildcards.seed,
            template=wildcards.template,
            lab=sorted(
                glob_wildcards(
                    "results/diffparc/sub-" + wildcards.subject +
                    "/probtrackx_output/seeds_to_sub-" + wildcards.subject + "_label-{lab}_mask.nii.gz"
                ).lab
            )
        ),
    output:
        done = "results/diffparc/sub-{subject}/probtrackx_output_seed-{seed}_space-{template}_desc-Linear/DONE.txt"
    envmodules: 'ants'
    container: config['singularity']['neuroglia']
    threads: 32
    resources:
        mem_mb = 128000
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{template}_Linear.log'
    group: 'participant1'
    shell:
        r'''
        mkdir -p $(dirname {output.done})
        ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 \
        parallel --jobs {threads} \
            antsApplyTransforms -d 3 --interpolation Linear \
            -i {{1}} -o {{2}} \
            -r {input.ref} -t {input.warp} \
            &> {log} \
            ::: {params.in_files} \
            :::+ {params.out_files}

        echo "done" > {output.done}
        '''



# BF マスクをコピー
rule get_binary_template_seed:
    input: 
        seed = lambda wildcards: config['template_binary_mask'][wildcards.seed]
    output: 
        mask = bids(root='results/diffparc',template='{template}',label='{seed}',suffix='mask.nii.gz')
    group: 'group0'
    shell:
        'cp {input} {output}'
        
# seed を3×3×3 voxel で膨張
rule dilate_seed:
    input: 
        mask = bids(root='results/diffparc',template='{template}',label='{seed}',suffix='mask.nii.gz')
    output: 
        seed = bids(root='results/diffparc',template='{template}',label='{seed}',desc='dilatedsmoothed',suffix='mask.nii.gz'),
    container: config['singularity']['neuroglia']
    group: 'group0'
    shell:
        'c3d {input} -dilate 1 3x3x3vox -o {output}'
 
# 拡張済みseed を被験者空間へ変換
rule transform_to_subject:
    input:
        seed = bids(root='results/diffparc',template='{template}',label='{seed}',desc='dilatedsmoothed',suffix='mask.nii.gz'),
        invwarp =  config['ants_invwarp_nii'],
        ref = '/home/m-ehara/project/PNI/derivatives/sub-{subject}/ses-01/parc/sub-{subject}_ses-01_space-nativepro_T1w_atlas-schaefer-400.nii.gz'
    output:
        seed = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='mask.nii.gz'),
    envmodules: 'ants'
    container: config['singularity']['neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}.log'
    group: 'participant1'
    threads: 8
    shell:
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.seed} -o {output.seed} -r {input.ref}  -t {input.invwarp} &> {log}'

# dwi brain mask を1.6 mm にそろえる
rule resample_brainmask:
    input:
        dwi = "/home/m-ehara/project/PNI/derivatives/sub-{subject}/ses-01/dwi/sub-{subject}_ses-01_space-dwi_desc-brain_mask.nii.gz"
    params:
        seed_resolution = 1.6
    output:
        mask = "results/diffparc/sub-{subject}/sub-{subject}_label-brain_mask.nii.gz",
        mask_res = "results/diffparc/sub-{subject}/sub-{subject}_label-brain_res-dwi_mask.nii.gz",
    log: "logs/resample_brainmask/sub-{subject}.log"
    shell:
        "fslmaths {input.dwi} -bin {output.mask} && "
        "mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}"

# targets atlas（皮質Schaefer400 ）の解像度を1.6 mm へリサンプル
rule resample_targets:
    input:
        mask_res = "results/diffparc/sub-{subject}/sub-{subject}_label-brain_res-dwi_mask.nii.gz",
        targets  = "/home/m-ehara/project/PNI/derivatives/sub-{subject}/ses-01/parc/sub-{subject}_ses-01_space-nativepro_T1w_atlas-schaefer-400.nii.gz"
    output:
        targets_res = "results/diffparc/sub-{subject}/sub-{subject}_space-individual_label-schaefer-400_res-dwi_dseg.nii.gz"
    container: config['singularity']['neuroglia']
    log: "logs/resample_targets/sub-{subject}.log"
    group: "participant1"
    shell:
        "reg_resample -ref {input.mask_res} -flo {input.targets} -res {output.targets_res} -inter 0 &> {log}"

# 被験者seed をdwi 解像度へリサンプル
rule resample_seed:
    input:
        mask_res = "results/diffparc/sub-{subject}/sub-{subject}_label-brain_res-dwi_mask.nii.gz",
        seed = "results/diffparc/sub-{subject}/sub-{subject}_space-individual_label-fullBF_from-MNI152NLin6Asym_mask.nii.gz"
    output:
        seed_res = "results/diffparc/sub-{subject}/sub-{subject}_space-individual_label-fullBF_from-MNI152NLin6Asym_res-dwi_mask.nii.gz"
    log:
        "logs/resample_seed/MNI152NLin6Asym_sub-{subject}_fullBF.log"
    group: "participant1"
    shell:
        "reg_resample -ref {input.mask_res} "
        "-flo {input.seed} "
        "-res {output.seed_res} "
        "-inter 0 &> {log}"

# target atlas の中で存在するラベル番号を400 個探して切り出す
rule split_targets:
    input:
        targets = "results/diffparc/sub-{subject}/sub-{subject}_space-individual_label-schaefer-400_res-dwi_dseg.nii.gz"
    output:
        target_seg_dir = directory("results/diffparc/sub-{subject}/targets")
    log:
        "logs/split_targets/sub-{subject}.log"
    shell:
        r"""
        mkdir -p {output.target_seg_dir}

        # 実在ラベル取得（背景0除外）
        labels=$(fslstats {input.targets} -H 3000 0 3000 \
                 | awk '$1>0 && NR-1!=0 {{print NR-1}}')

        for l in $labels; do
            fslmaths {input.targets} -thr $l -uthr $l -bin \
              {output.target_seg_dir}/sub-{wildcards.subject}_label-${{l}}_mask.nii.gz
        done
        """

# sprit_targetsで得られたパスを書く
import glob

rule gen_targets_txt:
    input:
        target_seg_dir = "results/diffparc/sub-{subject}/targets"
    output:
        target_txt = "results/diffparc/sub-{subject}/targets.txt"
    log:
        "logs/get_targets_txt/sub-{subject}.log"
    run:
        # targets フォルダ内の実ファイルを拾う
        files = sorted(glob.glob(f"{input.target_seg_dir}/sub-{wildcards.subject}_label-*_mask.nii.gz"))

        if len(files) == 0:
            raise ValueError(f"No target masks found in {input.target_seg_dir}")

        # 絶対パスで書きたいならここで変換
        files = ["/home/m-ehara/project/diffparc-smk/" + f for f in files]

        with open(output.target_txt, "w") as f:
            for path in files:
                f.write(path + "\n")


# BF から皮質のストリームライン到達を出す
rule run_probtrack:
    input:
        seed_res   = "results/diffparc/sub-{subject}/sub-{subject}_space-individual_label-fullBF_from-MNI152NLin6Asym_res-dwi_mask.nii.gz",
        target_txt = "results/diffparc/sub-{subject}/targets.txt",
        mask       = "results/diffparc/sub-{subject}/sub-{subject}_label-brain_res-dwi_mask.nii.gz",
        target_seg_dir = "results/diffparc/sub-{subject}/targets"
    params:
        bedpost_merged = join(config['fsl_bedpost_dir'], config['bedpost_merged_prefix']),
        probtrack_opts = config['probtrack']['opts'],
        nsamples       = config['probtrack']['nsamples']
    output:
        probtrack_dir = directory("results/diffparc/sub-{subject}/sub-{subject}_label-fullBF_from-MNI152NLin6Asym_probtrack")
    log:
        "logs/run_probtrack/sub-{subject}_fullBF_from-MNI152NLin6Asym.log"
    threads: 32
    shell:
        """
        mkdir -p {output.probtrack_dir}
        probtrackx2 \
            --samples={params.bedpost_merged} \
            --mask={input.mask} \
            --seed={input.seed_res} \
            --targetmasks={input.target_txt} \
            --seedref={input.seed_res} \
            --nsamples={params.nsamples} \
            --dir={output.probtrack_dir} \
            {params.probtrack_opts} -V 2 &> {log}
        """


