import argparse
import numpy as np
import nibabel as nib
from dipy.core.gradients import gradient_table
from dipy.align import register_dwi_to_template
import ants
from tqdm import tqdm
import shutil

def run(dwi_path, atlas_path, warp_file=None, affine_file=None, rev_warp_file=None, rev_affine_file=None):
    """Register a diffusion-weighted image to a reference atlas using ANTs.
    
    This function performs registration between a diffusion-weighted image (DWI) and 
    a reference atlas using ANTs' SyNRA transform, which includes both linear (rigid + affine) 
    and nonlinear (SyN) components. The transformation files can optionally be saved
    for later use.
    
    Args:
        dwi_path (str): Path to the input DWI image (NIfTI file).
        atlas_path (str): Path to the reference atlas image (NIfTI file).
        warp_file (str, optional): Path to save the forward warp field. Defaults to None.
        affine_file (str, optional): Path to save the forward affine transform. Defaults to None.
        rev_warp_file (str, optional): Path to save the reverse warp field. Defaults to None.
        rev_affine_file (str, optional): Path to save the reverse affine transform. Defaults to None.
        
    Returns:
        None: The function saves transformation files to the specified paths but does not
        return any values.
        
    Notes:
        The function uses ANTsPy to perform the registration and assumes the DWI image
        is already preprocessed (bias-corrected, possibly denoised). The function saves
        both forward and inverse transformations if paths are provided.
    """
    # Read the main DWI file using ANTs
    dwi_ants = ants.image_read(dwi_path)
    dwi_data = dwi_ants.numpy()

    atlas_ants = ants.image_read(atlas_path)    
    atlas_data = atlas_ants.numpy()
    

    b0_ants = ants.from_numpy(
        dwi_data,
        origin=dwi_ants.origin[:3],
        spacing=dwi_ants.spacing[:3]
    )

    
    # 'SyN' transform includes both linear and nonlinear registration.
    transforms = ants.registration(fixed=atlas_ants, moving=b0_ants, type_of_transform="SyNRA")

    # The result of the registration is a dictionary containing, among other keys:
    # 'warpedmovout' and 'fwdtransforms' (list of transform paths generated).
    registered = transforms["warpedmovout"]

    # If specified, save the transform files
    # Typically, transforms["fwdtransforms"][0] is the warp field, and [1] is the affine.
    if warp_file:
        shutil.copyfile(transforms["fwdtransforms"][0], warp_file)
        print(f"Saved warp field as {warp_file}")
    if affine_file:
        shutil.copyfile(transforms["fwdtransforms"][1], affine_file)
        print(f"Saved affine transform as {affine_file}")
    if rev_warp_file:
        shutil.copyfile(transforms["invtransforms"][0], rev_warp_file)
        print(f"Saved reverse warp field as {rev_warp_file}")
    if rev_affine_file:
        shutil.copyfile(transforms["invtransforms"][1], rev_affine_file)
        print(f"Saved reverse affine transform as {rev_affine_file}")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform linear registration of bias-corrected DWI to an atlas."
    )
    parser.add_argument("--moving", type=str, required=True,
                        help="Path to the bias-corrected DWI image (NIfTI file).")
    parser.add_argument("--fixed", type=str, required=True,
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
    run(args.moving, args.fixed, args.warpfield, args.affine, args.rev_warpfield, args.rev_affine)
    print("Registration complete.")
    