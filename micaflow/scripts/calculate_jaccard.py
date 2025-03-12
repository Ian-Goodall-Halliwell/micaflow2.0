import csv
import nibabel as nib
import argparse
import numpy as np
import sys

def print_help_message():
    help_text = """
    ╔════════════════════════════════════════════════════════════════╗
    ║                    JACCARD INDEX CALCULATOR                    ║
    ╚════════════════════════════════════════════════════════════════╝
    
    This script calculates the Jaccard similarity index (intersection over union)
    between two segmentation volumes, either globally or for each ROI.
    
    REQUIRED ARGUMENTS:
      --input, -i     : Path to the first input segmentation volume (.nii.gz)
      --reference, -r : Path to the reference segmentation volume (.nii.gz)
      --output, -o    : Output path for the CSV file with Jaccard indices
    
    OPTIONAL ARGUMENTS:
      --mask, -m      : Optional mask to restrict comparison to a specific region
      --threshold, -t : Threshold value for probabilistic segmentations (default: 0.5)
    
    EXAMPLE USAGE:
      micaflow calculate_jaccard \\
        --input segmentation1.nii.gz \\
        --reference ground_truth.nii.gz \\
        --output jaccard_metrics.csv
      
      # With mask and custom threshold:
      micaflow calculate_jaccard \\
        --input segmentation1.nii.gz \\
        --reference ground_truth.nii.gz \\
        --output jaccard_metrics.csv \\
        --mask brain_mask.nii.gz \\
        --threshold 0.75
    
    NOTES:
    - For multi-label segmentations, the Jaccard index is computed for each label
    - Values range from 0 (no overlap) to 1 (perfect overlap)
    - A global Jaccard index is calculated across all labels
    
    """
    print(help_text)


def Overlap(volume1_path, volume2_path, mask_path=None):
    """
    Calculate Jaccard index between two segmented volumes.
    
    Args:
        volume1_path (str): Path to first volume
        volume2_path (str): Path to second volume
        mask_path (str, optional): Path to mask volume
        
    Returns:
        dict: Dictionary containing ROI-wise Jaccard indices
    """
    import numpy as np
    
    # Load volumes
    vol1_img = nib.load(volume1_path)
    vol2_img = nib.load(volume2_path)
    
    vol1_data = vol1_img.get_fdata()
    vol2_data = vol2_img.get_fdata()
    
    # Apply mask if provided
    if mask_path:
        mask_img = nib.load(mask_path)
        mask_data = mask_img.get_fdata().astype(bool)
        vol1_data = np.logical_and(vol1_data, mask_data)
        vol2_data = np.logical_and(vol2_data, mask_data)
    
    # Get unique ROIs (assuming ROIs are represented by integer values)
    roi_values = np.unique(np.where(vol1_data > 0, vol1_data, 0))
    roi_values = roi_values[roi_values > 0]  # Remove background (0)
    
    # Calculate Jaccard for each ROI
    roi_ji = []
    
    if len(roi_values) > 0:
        for roi in roi_values:
            roi1 = vol1_data == roi
            roi2 = vol2_data == roi
            
            intersection = np.logical_and(roi1, roi2).sum()
            union = np.logical_or(roi1, roi2).sum()
            
            # Calculate Jaccard index
            jaccard = intersection / union if union > 0 else 0
            roi_ji.append(jaccard)
        # Calculate global Jaccard if no ROIs found
        intersection = np.logical_and(vol1_data, vol2_data).sum()
        union = np.logical_or(vol1_data, vol2_data).sum()
        jaccard = intersection / union if union > 0 else 0
        roi_ji.append(jaccard)
    else:
        # Calculate global Jaccard if no ROIs found
        intersection = np.logical_and(vol1_data, vol2_data).sum()
        union = np.logical_or(vol1_data, vol2_data).sum()
        jaccard = intersection / union if union > 0 else 0
        roi_ji.append(jaccard)
    
    # Create a results object similar to nipype's Overlap
    class Results:
        class Outputs:
            def __init__(self, roi_ji):
                self.roi_ji = roi_ji
        
        def __init__(self, roi_ji):
            self.outputs = self.Outputs(roi_ji)
            
    return Results(roi_ji)

def main(image, reference, output_file, threshold=0.5, mask_path=None):
    # Apply threshold and use the new file paths

    # Use our custom Overlap function instead of nipype
    if mask_path:
        res = Overlap(image, reference, mask_path)
    else:
        res = Overlap(image, reference)

    # Print the number of ROIs
    num_rois = len(res.outputs.roi_ji)
    print("Number of ROIs:", num_rois)


    with open(output_file, "w", newline="") as file:
        csvwriter = csv.writer(file)
        csvwriter.writerow(["ROI", "Jaccard Index"])
        for i, ji in enumerate(res.outputs.roi_ji):
            csvwriter.writerow([i + 1, ji])


if __name__ == "__main__":
    # Check if no arguments were provided or help was requested
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print_help_message()
        sys.exit(0)
    
    parser = argparse.ArgumentParser(description="Calculate overlap metrics between two volumes")
    parser.add_argument("--input", "-i", required=True, help="First input volume")
    parser.add_argument("--reference", "-r", required=True, help="Reference volume to compare against")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path")
    parser.add_argument("--mask", "-m", help="Optional mask volume")
    parser.add_argument("--threshold", "-t", type=float, default=0.5, help="Threshold value (default: 0.5)")
    
    args = parser.parse_args()
    
    main(args.input, args.reference, args.output, threshold=args.threshold, mask_path=args.mask)