import argparse
import sys
from dipy.reconst.dti import TensorModel
from dipy.core.gradients import gradient_table
import nibabel as nib

def print_help_message():
    help_text = """
    ╔════════════════════════════════════════════════════════════════╗
    ║                DIFFUSION TENSOR METRICS (FA/MD)                ║
    ╚════════════════════════════════════════════════════════════════╝
    
    This script computes Fractional Anisotropy (FA) and Mean Diffusivity (MD)
    maps from diffusion-weighted images using the diffusion tensor model.
    
    REQUIRED ARGUMENTS:
      --input      : Path to the input DWI image (.nii.gz)
      --mask       : Path to the brain mask image (.nii.gz)
      --bval       : Path to the b-values file (.bval)
      --bvec       : Path to the b-vectors file (.bvec)
      --output-fa  : Output path for the FA map (.nii.gz)
      --output-md  : Output path for the MD map (.nii.gz)
    
    EXAMPLE USAGE:
      micaflow compute_fa_md \\
        --input corrected_dwi.nii.gz \\
        --mask brain_mask.nii.gz \\
        --bval dwi.bval \\
        --bvec dwi.bvec \\
        --output-fa fa.nii.gz \\
        --output-md md.nii.gz
    
    NOTES:
    - FA (Fractional Anisotropy) values range from 0 (isotropic) to 1 (anisotropic)
    - MD (Mean Diffusivity) measures the overall magnitude of diffusion
    - Processing requires a brain mask to exclude non-brain regions
    
    """
    print(help_text)

# ----- Function: FA/MD Estimation -----
def compute_fa_md(bias_corr_path, mask_path, moving_bval, moving_bvec, fa_path, md_path):
    """Compute Fractional Anisotropy (FA) and Mean Diffusivity (MD) maps from diffusion-weighted images.
    
    This function takes a bias-corrected diffusion-weighted image (DWI) and a brain mask,
    creates a diffusion tensor model, and calculates FA and MD maps. The resulting
    maps are saved as NIfTI files at the specified output paths.
    
    Args:
        bias_corr_path (str): Path to the bias-corrected DWI image (NIfTI file).
        mask_path (str): Path to the brain mask image (NIfTI file).
        moving_bval (str): Path to the b-values file (.bval).
        moving_bvec (str): Path to the b-vectors file (.bvec).
        fa_path (str): Output path for the fractional anisotropy (FA) map.
        md_path (str): Output path for the mean diffusivity (MD) map.
        
    Returns:
        tuple: A tuple containing two strings (fa_path, md_path) - the paths to the 
              saved FA and MD NIfTI files.
    """
    bias_corr = nib.load(bias_corr_path)
    mask = nib.load(mask_path)
    masked_data = bias_corr.get_fdata() * mask.get_fdata()[..., None]
    gtab = gradient_table(moving_bval, moving_bvec)
    tensor_model = TensorModel(gtab)
    tensor_fit = tensor_model.fit(masked_data)
    fa = tensor_fit.fa
    md = tensor_fit.md
    nib.save(nib.Nifti1Image(fa, bias_corr.affine), fa_path)
    nib.save(nib.Nifti1Image(md, bias_corr.affine), md_path)
    return fa_path, md_path

if __name__ == "__main__":
    # Check if no arguments were provided or help was requested
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print_help_message()
        sys.exit(0)
    
    parser = argparse.ArgumentParser(
        description="Compute FA and MD maps using bias-corrected DWI and a brain mask."
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Path to the bias-corrected DWI image (NIfTI file).")
    parser.add_argument("--mask", type=str, required=True,
                        help="Path to the brain mask image (NIfTI file).")
    parser.add_argument("--bval", type=str, required=True,
                        help="Path to the bvals file.")
    parser.add_argument("--bvec", type=str, required=True,
                        help="Path to the bvecs file.")
    parser.add_argument("--output-fa", type=str, required=True,
                        help="Output path for the FA map.")
    parser.add_argument("--output-md", type=str, required=True,
                        help="Output path for the MD map.")
    args = parser.parse_args()
    
    fa_path, md_path = compute_fa_md(args.input, args.mask, args.bval, args.bvec, args.output_fa, args.output_md)
    print("FA map saved as:", fa_path)
    print("MD map saved as:", md_path)