import subprocess
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform N4 Bias Field Correction")
    parser.add_argument("--input", "-i", required=True, help="Input image file")
    parser.add_argument(
        "--output", "-o", required=True, help="Output corrected image file"
    )
    args = parser.parse_args()
    input_abs_path = os.path.abspath(args.input)

    subprocess.run(
        "python3 scripts/HD_BET/entry_point.py -i "
        + input_abs_path
        + " -o "
        + args.output
        + " --save_bet_mask",  # + " -device cpu",
        shell=True,
    )
