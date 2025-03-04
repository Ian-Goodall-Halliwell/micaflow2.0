import argparse
import numpy as np
import nibabel as nib
from dipy.core.gradients import gradient_table
from dipy.align import register_dwi_to_template
import ants
from tqdm import tqdm
import shutil

def run(dwi_path, atlas_path, warp_file=None, affine_file=None, rev_warp_file=None, rev_affine_file=None):
    """
    Replace the previous motion correction logic with the QuickSyN-based 
    registration you provided for each volume in the DWI.
    """
    # Read the main DWI file using ANTs
    dwi_ants = ants.image_read(dwi_path)
    dwi_data = dwi_ants.numpy()

    atlas_ants = ants.image_read(atlas_path)    
    atlas_data = atlas_ants.numpy()
    
    # B0 is assumed to be the first volume (index 0)
    b0_data = dwi_data[..., 0]
    b0_ants = ants.from_numpy(
        b0_data,
        origin=dwi_ants.origin[:3],
        spacing=dwi_ants.spacing[:3]
    )

 
    # Rigid registration
    rigid_reg = ants.registration(
        fixed=atlas_ants,
        moving=b0_ants,
        type_of_transform='QuickRigid'
    )
    # Non-linear registration (SyNOnly) using the rigid transform as initial
    quicksyn_reg = ants.registration(
        fixed=atlas_ants,
        moving=rigid_reg['warpedmovout'],
        initial_transform=rigid_reg['fwdtransforms'][0],
        type_of_transform='SyNOnly'
    )

  
    # If specified, save the transform files
    # Typically, transforms["fwdtransforms"][0] is the warp field, and [1] is the affine.
    if warp_file:
        shutil.copyfile(quicksyn_reg["fwdtransforms"][0], warp_file)
        print(f"Saved warp field as {warp_file}")
    if affine_file:
        shutil.copyfile(rigid_reg["fwdtransforms"][0], affine_file)
        print(f"Saved affine transform as {affine_file}")
    if rev_warp_file:
        shutil.copyfile(quicksyn_reg["invtransforms"][0], rev_warp_file)
        print(f"Saved reverse warp field as {rev_warp_file}")
    if rev_affine_file:
        shutil.copyfile(rigid_reg["invtransforms"][0], rev_affine_file)
        print(f"Saved reverse affine transform as {rev_affine_file}")


# ----- Function: Linear Registration -----
def run_linear_registration(bias_corr_path, moving_bval, moving_bvec, atlas, affine_path):
    # Linear registration to atlas
    pipeline = ["center_of_mass", "translation", "rigid", "affine"]
    level_iters = [500, 100, 50]    # Adjusted parameters
    sigmas = [4.0, 2.0, 1.0]
    factors = [8, 4, 2]
    
    bias_corr = nib.load(bias_corr_path)
    MNI_atlas = nib.load(atlas)
    xformed_dwi, reg_affine = register_dwi_to_template(
        dwi=bias_corr,
        gtab=gradient_table(moving_bval, moving_bvec),
        template=MNI_atlas,
        reg_method="aff",
        nbins=32,
        metric='MI',
        pipeline=pipeline,
        level_iters=level_iters,
        sigmas=sigmas,
        factors=factors)
    # Save the affine matrix for later use
    np.savetxt(affine_path, reg_affine)
    fixed_xformed_path = "registered.nii.gz"
    fixed_xformed_dwi = nib.Nifti1Image(xformed_dwi, MNI_atlas.affine)
    nib.save(fixed_xformed_dwi, fixed_xformed_path)
    return fixed_xformed_path, reg_affine, MNI_atlas

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform linear registration of bias-corrected DWI to an atlas."
    )
    parser.add_argument("--bias_corr", type=str, required=True,
                        help="Path to the bias-corrected DWI image (NIfTI file).")
    parser.add_argument("--atlas", type=str, required=True,
                        help="Path to the atlas image (NIfTI file).")
    parser.add_argument("--affine", type=str, required=True,
                        help="Path for the affine output.")
    parser.add_argument("--rev_affine", type=str, required=True,
                        help="Path for the affine output.")
    parser.add_argument("--warpfield", type=str, required=True,
                        help="Path for the affine output.")
    parser.add_argument("--rev_warpfield", type=str, required=True,
                        help="Path for the affine output.")
    
    args = parser.parse_args()
    run(args.bias_corr, args.atlas, args.warpfield, args.affine, args.rev_warpfield, args.rev_affine)
    print("Registration complete.")
    