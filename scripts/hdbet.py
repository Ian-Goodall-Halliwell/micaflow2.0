import subprocess
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform N4 Bias Field Correction")
    parser.add_argument("--input", "-i", required=True, help="Input image file")
    parser.add_argument(
        "--output", "-o", required=True, help="Output corrected image file"
    )
    args = parser.parse_args()
    subprocess.run(
        "hd-bet -i "
        + args.input
        + " -o "
        + args.output
        + " --save_bet_mask",  # + " -device cpu",
        shell=True,
    )
