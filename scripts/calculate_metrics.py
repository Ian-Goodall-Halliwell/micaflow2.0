import sys
import csv
import os
import shutil
import tempfile
from nipype.algorithms.metrics import Overlap
import nibabel as nib

def apply_threshold(image_path, threshold=0.5):
    img = nib.load(image_path)
    data = img.get_fdata()
    data[data < threshold] = 0
    # Write to a new file to avoid in-place modifications that may cause issues
    new_image_path = image_path.replace(".nii", "_thr.nii")
    nib.save(nib.Nifti1Image(data, img.affine), new_image_path)
    return new_image_path

def main(image, reference, output_file, threshold=0.5, mask_path=None):
    # Create a temporary directory
    temp_dir = tempfile.mkdtemp()

    try:
        # Copy necessary files into the temporary directory
        temp_image = os.path.join(temp_dir, os.path.basename(image))
        temp_reference = os.path.join(temp_dir, os.path.basename(reference))
        shutil.copy(image, temp_image)
        shutil.copy(reference, temp_reference)

        if mask_path:
            temp_mask = os.path.join(temp_dir, os.path.basename(mask_path))
            shutil.copy(mask_path, temp_mask)

        # Change the current working directory to the temporary directory
        original_cwd = os.getcwd()
        os.chdir(temp_dir)

        # Apply threshold and use the new file paths
        image_thr = apply_threshold(temp_image, threshold)
        reference_thr = apply_threshold(temp_reference, threshold)

        overlap = Overlap()
        overlap.inputs.volume1 = image_thr
        overlap.inputs.volume2 = reference_thr

        if mask_path:
            # Optionally process mask similarly if thresholding is required
            mask_thr = apply_threshold(temp_mask, threshold)
            overlap.inputs.mask_volume = mask_thr

        res = overlap.run()

        # Print the number of ROIs
        num_rois = len(res.outputs.roi_ji)
        print("Number of ROIs:", num_rois)

        temp_output_file = os.path.join(temp_dir, "output.csv")
        with open(temp_output_file, "w", newline="") as file:
            csvwriter = csv.writer(file)
            csvwriter.writerow(["ROI", "Jaccard Index"])
            for i, ji in enumerate(res.outputs.roi_ji):
                csvwriter.writerow([i + 1, ji])

        # Copy the output CSV file to the desired location
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        shutil.copy(temp_output_file, output_file)
    finally:
        # Change back to the original working directory
        os.chdir(original_cwd)
        # Remove the temporary directory and its contents
        shutil.rmtree(temp_dir)

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(
            "Usage: python calculate_metrics.py <volume1> <volume2> <output_csv> [<mask_path>]"
        )
        sys.exit(1)

    volume1 = sys.argv[1]
    volume2 = sys.argv[2]
    output_csv = sys.argv[3]
    mask_path = sys.argv[4] if len(sys.argv) == 5 else None

    main(volume1, volume2, output_csv, mask_path=mask_path)