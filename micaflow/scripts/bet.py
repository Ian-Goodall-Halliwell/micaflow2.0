"""HD-BET (High-Definition Brain Extraction Tool) wrapper script.

This script provides a simplified command-line interface to the HD-BET brain extraction
tool, which performs accurate skull stripping on brain MR images using a deep learning approach.
It supports both CPU and GPU execution modes.

The script is a wrapper around the HD-BET entry_point.py script that simplifies the interface
and handles path resolution.

Example:
    python hdbet.py --input t1w.nii.gz --output t1w_brain.nii.gz
    python hdbet.py --input t1w.nii.gz --output t1w_brain.nii.gz --cpu

"""

import subprocess
import argparse
import os
import shutil
import sys

def print_help_message():
    help_text = """
    ╔════════════════════════════════════════════════════════════════╗
    ║                           HD-BET                               ║
    ╚════════════════════════════════════════════════════════════════╝
    
    This script performs brain extraction (skull stripping) on MRI images 
    using the HD-BET deep learning tool. It accurately segments the brain 
    from surrounding tissues.
    
    REQUIRED ARGUMENTS:
      --input, -i      : Path to the input MR image (.nii.gz)
      --output, -o     : Path for the output brain-extracted image (.nii.gz)
      --output-mask, -m: Path for the output brain mask (.nii.gz)
    
    OPTIONAL ARGUMENTS:
      --cpu            : Use CPU instead of GPU for computation (slower but works without CUDA)
    
    EXAMPLE USAGE:
      micaflow bet \\
        --input t1w.nii.gz \\
        --output t1w_brain.nii.gz \\
        --output-mask t1w_brain_mask.nii.gz

      # Run on CPU instead of GPU:
      micaflow bet \\
        --input t1w.nii.gz \\
        --output t1w_brain.nii.gz \\
        --output-mask t1w_brain_mask.nii.gz \\
        --cpu
    
    NOTES:
    - GPU acceleration is used by default for faster processing
    - The output is a brain-extracted image and a binary brain mask
    - HD-BET uses a 3D U-Net architecture for robust brain extraction
    
    """
    print(help_text)


if __name__ == "__main__":
    # Check if no arguments were provided or help was requested
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print_help_message()
        sys.exit(0)
        
    parser = argparse.ArgumentParser(description="Perform brain extraction using HD-BET")
    parser.add_argument("--input", "-i", required=True, help="Input MR image file")
    parser.add_argument(
        "--output", "-o", required=True, help="Output brain-extracted image file"
    )
    parser.add_argument(
        "--output-mask", "-m", required=True, help="Output brain-extracted mask image file"
    )
    parser.add_argument("--cpu", action="store_true", help="Use CPU instead of GPU")
    args = parser.parse_args()
    input_abs_path = os.path.abspath(args.input)

    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Create path to HD_BET entry_point.py relative to this script's location
    hdbet_script = os.path.join(script_dir, "HD_BET", "entry_point.py")

    subprocess.run(
        f"python3 {hdbet_script} -i {input_abs_path} -o {args.output} --save_bet_mask{' -device cpu --disable_tta' if args.cpu else ''}",
        shell=True,
    )
    shutil.move(args.output.replace(".nii.gz", "") + "_bet" + ".nii.gz", args.output_mask)