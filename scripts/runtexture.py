from texturepipeline import noelTexturesPy
import argparse


def run_texture_pipeline(input, mask, output_dir):
    """Run the neuroimaging texture feature extraction pipeline.
    
    This function initializes and executes a texture analysis pipeline on a neuroimaging volume.
    The pipeline computes various texture features (e.g., gradient magnitude, relative intensity,
    local binary patterns) from the input image within the regions defined by the mask.
    Results are saved to the specified output directory.
    
    Parameters
    ----------
    input : str
        Path to the input image file (typically a preprocessed MRI volume).
    mask : str
        Path to the binary mask file that defines regions of interest for texture analysis.
    output_dir : str
        Directory where the computed texture feature maps will be saved.
    
    Returns
    -------
    None
        The function saves texture feature maps to the output directory but does not return values.
        
    Notes
    -----
    The function relies on the noelTexturesPy class which implements multiple texture
    feature extraction algorithms specifically designed for neuroimaging data.
    """
    pipeline = noelTexturesPy(
        id='textures',
        output_dir=output_dir,
        input=input,
        mask=mask,
    )
    pipeline.file_processor()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform N4 Bias Field Correction")
    parser.add_argument("--input", "-i", required=True, help="Input image file")
    parser.add_argument("--mask", "-m", required=True, help="Input mask file")
    parser.add_argument(
        "--output", "-o", required=True, help="Output corrected image file"
    )
    args = parser.parse_args()
    run_texture_pipeline(args.input, args.mask, args.output)
