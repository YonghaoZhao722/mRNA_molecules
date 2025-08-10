import os
import tifffile
import numpy as np
import json
from scipy.ndimage import binary_dilation
from skimage.morphology import disk

# Configuration parameters
add_noise = True  # Set to False to disable Gaussian noise addition

def create_background_noise(image_slice, mask_single, band_width=5, mask=None):
    """Create background noise based on neighboring regions"""
    dilated_mask = binary_dilation(mask_single, structure=disk(band_width))
    band_region = dilated_mask & ~mask_single
    if np.sum(band_region) > 0:
        band_values = image_slice[band_region]
        mean_val = np.mean(band_values)
        std_val = np.std(band_values)
        noise_std = 0.2 * std_val
        return mean_val, noise_std
    else:
        if mask is not None:
            background_region = mask == 0
            if np.sum(background_region) > 0:
                bg_values = image_slice[background_region]
                mean_val = np.mean(bg_values)
                std_val = np.std(bg_values)
                return mean_val, 0.2 * std_val
        return 0, 0

def extract_cells_single(mask_path, image_path, output_dir):
    """Extract individual cells from a single mask and image pair"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Load mask and image
    mask = tifffile.imread(mask_path)
    image = tifffile.imread(image_path)
    
    # Validate dimensions
    if mask.shape != image.shape[1:]:
        raise ValueError(f"Mask xy size does not match image: mask {mask.shape}, image {image.shape[1:]}")
    
    # Get unique labels (cell IDs) from mask
    labels = np.unique(mask)
    labels = labels[labels != 0]  # Remove background (0)
    
    print(f"Found {len(labels)} cells in mask")
    
    # Sort cells by y position (top to bottom)
    cell_info = []
    for label in labels:
        mask_single = (mask == label)
        coords = np.argwhere(mask_single)
        y_min, x_min = coords.min(axis=0)
        cell_info.append((label, y_min))
    
    cell_info.sort(key=lambda x: x[1])
    
    # Create label mapping (renumber cells 1, 2, 3, ...)
    label_mapping = {}
    for i, (original_label, y_min) in enumerate(cell_info):
        label_mapping[original_label] = i + 1
    
    # Extract each cell
    coordinate_mapping = {}
    for original_label, new_id in label_mapping.items():
        print(f"Processing cell {new_id} (original label {original_label})")
        
        # Get cell boundaries
        mask_single = (mask == original_label)
        coords = np.argwhere(mask_single)
        y_min, x_min = coords.min(axis=0)
        y_max, x_max = coords.max(axis=0)
        
        # Calculate crop region (square bounding box)
        h = y_max - y_min + 1
        w = x_max - x_min + 1
        side = max(h, w)
        
        y_center = (y_min + y_max) // 2
        x_center = (x_min + x_max) // 2
        half = side // 2
        
        # Ensure crop region is within image bounds
        y1 = max(0, y_center - half)
        y2 = min(mask.shape[0], y1 + side)
        x1 = max(0, x_center - half)
        x2 = min(mask.shape[1], x1 + side)
        
        # Adjust if needed to maintain square
        if y2 - y1 < side:
            y1 = max(0, y2 - side)
        if x2 - x1 < side:
            x1 = max(0, x2 - side)
        
        # Store coordinate mapping
        coordinate_mapping[f'cell_{new_id:03d}'] = {
            'cell_id': int(new_id),
            'original_mask_label': int(original_label),
            'mask_bbox': {
                'y_min': int(y_min),
                'y_max': int(y_max),
                'x_min': int(x_min),
                'x_max': int(x_max)
            },
            'crop_region': {
                'y_start': int(y1),
                'y_end': int(y2),
                'x_start': int(x1),
                'x_end': int(x2),
                'y_offset': int(y1),
                'x_offset': int(x1),
                'crop_size': int(y2 - y1)
            },
            'center': {
                'y_center': int(y_center),
                'x_center': int(x_center)
            }
        }
        
        # Crop mask and image
        mask_crop = mask[y1:y2, x1:x2]
        image_crop = image[:, y1:y2, x1:x2]
        mask_single_crop = mask_crop == original_label
        
        # Process each z-slice
        result = np.zeros_like(image_crop, dtype=np.uint16)
        for z in range(image.shape[0]):
            current_slice = image_crop[z]
            result_slice = np.zeros_like(current_slice, dtype=np.uint16)
            
            if add_noise:
                # Add Gaussian noise to background
                mean_noise, std_noise = create_background_noise(current_slice, mask_single_crop, mask=mask_crop)
                noise_shape = current_slice.shape
                noise = np.random.normal(mean_noise, std_noise, noise_shape)
                noise = np.clip(noise, 0, 65535).astype(np.uint16)
                result_slice[:] = noise
            else:
                # Keep background as zeros (no noise)
                result_slice[:] = 0
            
            # Copy cell pixels
            result_slice[mask_single_crop] = current_slice[mask_single_crop]
            result[z] = result_slice
        
        # Save extracted cell
        output_path = os.path.join(output_dir, f'cell_{new_id:03d}.tif')
        tifffile.imwrite(output_path, result)
        print(f"Saved: {output_path}")
    
    # Save coordinate mapping
    mapping_file = os.path.join(output_dir, 'coordinate_mapping.json')
    with open(mapping_file, 'w') as f:
        json.dump(coordinate_mapping, f, indent=2)
    print(f"Saved coordinate mapping: {mapping_file}")

if __name__ == '__main__':
    # Input files
    mask_path = '/Volumes/ExFAT/deconv/3_s1/3_s1_mask.tif'
    image_path = '/Volumes/ExFAT/deconv/3_s1/3_s1_RL30.tif'
    output_dir = '/Volumes/ExFAT/deconv/3_s1/extracted_cells_RL30'
    
    print(f"Gaussian noise addition: {'Enabled' if add_noise else 'Disabled'}")
    print(f"Mask file: {mask_path}")
    print(f"Image file: {image_path}")
    print(f"Output directory: {output_dir}")
    
    # Check if input files exist
    if not os.path.exists(mask_path):
        print(f"Error: Mask file not found: {mask_path}")
        exit(1)
    
    if not os.path.exists(image_path):
        print(f"Error: Image file not found: {image_path}")
        exit(1)
    
    try:
        extract_cells_single(mask_path, image_path, output_dir)
        print(f"Processing completed successfully!")
        print(f"Extracted cells saved to: {output_dir}")
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()