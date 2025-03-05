import argparse
import nibabel as nib
from dipy.denoise.patch2self import patch2self
from dipy.io.gradients import read_bvals_bvecs


# ----- Function: Denoise -----
def run_denoise(moving, moving_bval, moving_bvec, output):
    moving_image = nib.load(moving)
    moving_bval_value, moving_bvec_value = read_bvals_bvecs(moving_bval, moving_bvec)
    denoised = patch2self(
        moving_image.get_fdata(),
        moving_bval_value,
        model="ols",
        shift_intensity=True,
        clip_negative_vals=False,
        b0_threshold=50,
        b0_denoising=False,
    )

    nib.save(nib.Nifti1Image(denoised, moving_image.affine), output)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Denoise a DWI image using patch2self."
    )
    parser.add_argument(
        "--moving",
        type=str,
        required=True,
        help="Path to the input DWI image (NIfTI file).",
    )
    parser.add_argument(
        "--bval", type=str, required=True, help="Path to the bvals file."
    )
    parser.add_argument(
        "--bvec", type=str, required=True, help="Path to the bvecs file."
    )
    parser.add_argument("--output", type=str, required=True, help="output path")

    args = parser.parse_args()
    output_path = run_denoise(args.moving, args.bval, args.bvec, args.output)
    print("Denoised image saved as:", output_path)
