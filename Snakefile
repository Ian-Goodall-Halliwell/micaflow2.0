import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Default parameters
SUBJECT = config.get("subject", "")
SESSION = config.get("session", "")
OUT_DIR = config.get("out_dir", "")
THREADS = config.get("threads", 1)
DATA_DIRECTORY = config.get("data_directory", "")
RUN_DWI = config.get("run_dwi", True)
CLEANUP = config.get("cleanup", True)

# Atlas paths
ATLAS_DIR = "atlas"
ATLAS = os.path.join(ATLAS_DIR, "mni_icbm152_t1_tal_nlin_sym_09a.nii")
ATLAS_MASK = os.path.join(ATLAS_DIR, "mni_icbm152_t1_tal_nlin_sym_09a_mask.nii")
ATLAS_SEG = os.path.join(ATLAS_DIR, "mni_icbm152_t1_tal_nlin_sym_09a_seg.nii")

def get_final_output():
    outputs = []
    # Add texture outputs
    for modality in ["T1w", "FLAIR"]:
        outputs.extend([
            f"{OUT_DIR}/{SUBJECT}/{SESSION}/textures/{SUBJECT}_{SESSION}_{modality}_textures_output_gradient_magnitude.nii",
            f"{OUT_DIR}/{SUBJECT}/{SESSION}/textures/{SUBJECT}_{SESSION}_{modality}_textures_output_relative_intensity.nii"
        ])
    # Add metrics outputs
    for modality in ["T1w", "FLAIR"]:
        outputs.append(f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/{modality}_{SUBJECT}_{SESSION}_jaccard.csv")
    
    if RUN_DWI:
        outputs.extend([
            f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/fa_registered.nii.gz",
            f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/md_registered.nii.gz"
        ])
    
    return outputs

rule all:
    input:
        get_final_output()

rule skull_strip:
    input:
        image = lambda wildcards: f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/anat/{SUBJECT}_{SESSION}_{'run-1_' if wildcards.modality == 'T1w' else ''}{wildcards.modality}.nii.gz"
    output:
        brain = f"{OUT_DIR}/{{modality}}_hdbet.nii.gz",
        mask = f"{OUT_DIR}/{{modality}}_hdbet_bet.nii.gz"
    conda:
        "envs/micaflow.yml"
    shell:
        "python3 scripts/hdbet.py --input {input.image} --output {output.brain}"

rule bias_field_correction:
    input:
        image = rules.skull_strip.output.brain,
        mask = rules.skull_strip.output.mask
    output:
        corrected = f"{OUT_DIR}/{SUBJECT}/{SESSION}/anat/{SUBJECT}_{SESSION}_desc-N4_{{modality}}.nii.gz"
    conda:
        "envs/micaflow.yml"
    shell:
        "python3 scripts/N4BiasFieldCorrection.py -i {input.image} -o {output.corrected} -m {input.mask}"

rule synthseg_t1w:
    input:
        image = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/anat/{SUBJECT}_{SESSION}_run-1_T1w.nii.gz"
    output:
        seg = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/T1w_parcellation.nii.gz"
    conda:
        "envs/micaflow.yml"
    threads: THREADS
    shell:
        "python3 scripts/run_synthseg.py --i {input.image} --o {output.seg} --parc --fast --threads {threads}"

rule synthseg_flair:
    input:
        image = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/anat/{SUBJECT}_{SESSION}_FLAIR.nii.gz",
        t1w_seg = rules.synthseg_t1w.output.seg
    output:
        seg = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/FLAIR_parcellation.nii.gz"
    conda:
        "envs/micaflow.yml"
    threads: THREADS
    shell:
        "python3 scripts/run_synthseg.py --i {input.image} --o {output.seg} --parc --fast --threads {threads}"

rule registration_t1w:
    input:
        fixed = rules.synthseg_t1w.output.seg,
        moving = rules.synthseg_flair.output.seg
    output:
        warped = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_FLAIR_space-T1w.nii.gz",
        fwd_field = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-FLAIR_to-T1w_fwdfield.nii.gz",
        bak_field = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-FLAIR_to-T1w_bakfield.nii.gz",
        fwd_affine = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-FLAIR_to-T1w_fwdaffine.mat",
        bak_affine = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-FLAIR_to-T1w_bakaffine.mat"
    conda:
        "envs/micaflow.yml"
    shell:
        """
        python3 scripts/coregister.py \
            --fixed-file {input.fixed} \
            --moving-file {input.moving} \
            --out-file {output.warped} \
            --warp-file {output.fwd_field} \
            --affine-file {output.fwd_affine} \
            --rev-warp-file {output.bak_field} \
            --rev-affine-file {output.bak_affine}
        """

rule registration_mni152:
    input:
        image = rules.synthseg_t1w.output.seg,
        fixed = ATLAS_SEG
    output:
        warped = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_T1w_space-MNI152.nii.gz",
        fwd_field = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-T1w_to-MNI152_fwdfield.nii.gz",
        bak_field = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-T1w_to-MNI152_bakfield.nii.gz",
        fwd_affine = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-T1w_to-MNI152_fwdaffine.mat",
        bak_affine = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-T1w_to-MNI152_bakaffine.mat"
    conda:
        "envs/micaflow.yml"
    shell:
        """
        python3 scripts/coregister.py \
            --fixed-file {input.fixed} \
            --moving-file {input.image} \
            --out-file {output.warped} \
            --warp-file {output.fwd_field} \
            --affine-file {output.fwd_affine} \
            --rev-warp-file {output.bak_field} \
            --rev-affine-file {output.bak_affine}
        """

rule apply_warp_flair_to_t1w:
    input:
        moving = rules.bias_field_correction.output.corrected,
        warp = rules.registration_t1w.output.fwd_field,
        affine = rules.registration_t1w.output.fwd_affine,
        reference = rules.bias_field_correction.output.corrected
    output:
        warped = f"{OUT_DIR}/{SUBJECT}/{SESSION}/anat/{SUBJECT}_{SESSION}_space-T1w_FLAIR.nii.gz"
    conda:
        "envs/micaflow.yml"
    shell:
        """
        python3 scripts/use_warp.py \
            --moving {input.moving} \
            --reference {input.reference} \
            --affine {input.affine} \
            --warp {input.warp} \
            --out {output.warped}
        """

rule apply_warp_to_mni:
    input:
        moving = lambda wildcards: rules.bias_field_correction.output.corrected if wildcards.modality == "T1w" else rules.apply_warp_flair_to_t1w.output.warped,
        warp = rules.registration_mni152.output.fwd_field,
        affine = rules.registration_mni152.output.fwd_affine,
        reference = ATLAS
    output:
        warped = f"{OUT_DIR}/{SUBJECT}/{SESSION}/anat/{SUBJECT}_{SESSION}_space-MNI152_{{modality}}.nii.gz"
    conda:
        "envs/micaflow.yml"
    shell:
        """
        python3 scripts/use_warp.py \
            --moving {input.moving} \
            --reference {input.reference} \
            --affine {input.affine} \
            --warp {input.warp} \
            --out {output.warped}
        """

rule run_texture:
    input:
        image = rules.apply_warp_to_mni.output.warped,
        mask = ATLAS_MASK
    output:
        gradient = f"{OUT_DIR}/{SUBJECT}/{SESSION}/textures/{SUBJECT}_{SESSION}_{{modality}}_textures_output_gradient_magnitude.nii",
        intensity = f"{OUT_DIR}/{SUBJECT}/{SESSION}/textures/{SUBJECT}_{SESSION}_{{modality}}_textures_output_relative_intensity.nii"
    conda:
        "envs/micaflow.yml"
    shell:
        """
        python3 scripts/runtexture.py \
            --input {input.image} \
            --mask {input.mask} \
            --output {OUT_DIR}/{SUBJECT}/{SESSION}/textures/{SUBJECT}_{SESSION}_{wildcards.modality}_textures_output
        """

rule calculate_metrics:
    input:
        image = rules.apply_warp_to_mni.output.warped,
        atlas = ATLAS
    output:
        metrics = f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/{{modality}}_{SUBJECT}_{SESSION}_jaccard.csv"
    conda:
        "envs/micaflow.yml"
    shell:
        """
        python3 scripts/calculate_metrics.py \
            {input.image} \
            {input.atlas} \
            {output.metrics}
        """

if RUN_DWI:
    rule dwi_denoise:
        input:
            moving = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.nii.gz",
            bval = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.bval",
            bvec = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.bvec"
        output:
            denoised = f"{OUT_DIR}/{SUBJECT}/{SESSION}/dwi/denoised_moving.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_denoise.py \
                --moving {input.moving} \
                --bval {input.bval} \
                --bvec {input.bvec}
            """

    rule dwi_motion_correction:
        input:
            denoised = rules.dwi_denoise.output.denoised,
            bval = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.bval",
            bvec = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.bvec"
        output:
            corrected = f"{OUT_DIR}/{SUBJECT}/{SESSION}/dwi/moving_motion_corrected.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_motioncorrection.py \
                --denoised {input.denoised} \
                --bval {input.bval} \
                --bvec {input.bvec}
            """

    rule dwi_topup:
        input:
            moving = rules.dwi_motion_correction.output.corrected,
            b0 = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_PA.nii.gz",
            b0_bval = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_PA.bval",
            b0_bvec = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_PA.bvec"
        output:
            warp = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/topup-warp-EstFieldMap.nii.gz",
            corrected = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/corrected_image.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/pyhysco.py \
                --data_image {input.moving} \
                --reverse_image {input.b0} \
                --output_name {output.corrected}
            """

    rule dwi_apply_topup:
        input:
            motion_corr = rules.dwi_motion_correction.output.corrected,
            warp = rules.dwi_topup.output.warp,
            affine = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.nii.gz"
        output:
            corrected = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/topup_corrected.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_applytopup.py \
                --motion_corr {input.motion_corr} \
                --warp {input.warp} \
                --affine {input.affine}
            """

    rule dwi_skull_strip:
        input:
            image = rules.dwi_topup.output.corrected
        output:
            mask = f"{OUT_DIR}/{SUBJECT}/{SESSION}/anat/DWI_hdbet_bet.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/hdbet.py \
                --input {input.image} \
                --output DWI_hdbet.nii.gz
            """

    rule dwi_bias_correction:
        input:
            image = rules.dwi_apply_topup.output.corrected,
            mask = rules.dwi_skull_strip.output.mask
        output:
            corrected = f"{OUT_DIR}/{SUBJECT}/{SESSION}/dwi/denoised_moving_corrected.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_biascorrection.py \
                --image {input.image} \
                --mask {input.mask}
            """

    rule synthseg_dwi:
        input:
            image = rules.dwi_topup.output.corrected,
            flair_seg = rules.synthseg_flair.output.seg
        output:
            seg = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/DWI_parcellation.nii.gz"
        conda:
            "envs/micaflow.yml"
        threads: THREADS
        shell:
            """
            python3 scripts/run_synthseg.py \
                --i {input.image} \
                --o {output.seg} \
                --parc \
                --fast \
                --threads {threads}
            """

    rule dwi_registration:
        input:
            moving = rules.synthseg_dwi.output.seg,
            fixed = rules.synthseg_t1w.output.seg
        output:
            fwd_field = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-DWI_to-T1w_fwdfield.nii.gz",
            fwd_affine = f"{OUT_DIR}/{SUBJECT}/{SESSION}/xfm/{SUBJECT}_{SESSION}_from-DWI_to-T1w_fwdaffine.mat"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_reg.py \
                --fixed {input.fixed} \
                --moving {input.moving} \
                --affine {output.fwd_affine} \
                --warpfield {output.fwd_field}
            """

    rule dwi_compute_fa_md:
        input:
            image = rules.dwi_bias_correction.output.corrected,
            mask = rules.dwi_topup.output.corrected,
            bval = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.bval",
            bvec = f"{DATA_DIRECTORY}/{SUBJECT}/{SESSION}/dwi/{SUBJECT}_{SESSION}_b700.bvec"
        output:
            fa = f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/fa_map.nii.gz",
            md = f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/md_map.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_compute_fa_md.py \
                --bias_corr {input.image} \
                --mask {input.mask} \
                --bval {input.bval} \
                --bvec {input.bvec}
            """

    rule dwi_fa_md_registration:
        input:
            fa = rules.dwi_compute_fa_md.output.fa,
            md = rules.dwi_compute_fa_md.output.md,
            atlas = rules.bias_field_correction.output.corrected,
            affine = rules.dwi_registration.output.fwd_affine,
            warp = rules.dwi_registration.output.fwd_field
        output:
            fa_reg = f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/fa_registered.nii.gz",
            md_reg = f"{OUT_DIR}/{SUBJECT}/{SESSION}/metrics/md_registered.nii.gz"
        conda:
            "envs/micaflow.yml"
        shell:
            """
            python3 scripts/dwi_fa_md_registration.py \
                --fa {input.fa} \
                --md {input.md} \
                --atlas {input.atlas} \
                --reg_affine {input.affine} \
                --mapping {input.warp} \
                --out_fa {output.fa_reg} \
                --out_md {output.md_reg}
            """ 