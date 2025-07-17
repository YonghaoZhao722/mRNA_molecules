import os
import tifffile
import numpy as np
import json
import csv
from skimage.measure import regionprops, label

# Configuration parameters
base_dir = r'Y333 ATP6 TIM50'
deconv_dir = os.path.join(base_dir, 'deconvolved_30')
mask_dir = os.path.join(base_dir, 'aligned_masks')
output_base = os.path.join(base_dir, 'extracted_cells_30_with_geometry')
os.makedirs(output_base, exist_ok=True)

# Microscope parameters
pixel_size_xy = 0.0645  # µm (adjust according to your microscope setup)
z_step = 0.2  # µm (confirm FITC z-stack layer spacing)

# Volume calculation modes
calculate_triaxial_ellipsoid = False  # Uses z-height estimation: V = (4/3)πabc
calculate_prolate_ellipsoid = True   # Rotation around major axis: V = (4/3)πa²b

def parse_deconv_filename(filename):
    name = os.path.splitext(filename)[0]
    if not name.startswith('deconv_'):
        return None, None
    name = name[len('deconv_'):]
    if '_w4FITC-100-' in name:
        prefix, sidx = name.split('_w4FITC-100-')
    else:
        parts = name.split('_s')
        if len(parts) < 2:
            return None, None
        prefix = '_s'.join(parts[:-1])
        sidx = 's' + parts[-1]
        return prefix, sidx
    if sidx.startswith('_'):
        sidx = sidx[1:]
    return prefix, sidx

def find_mask_file(prefix, sidx):
    mask_name = f'{prefix}_DIC_{sidx}.tif'
    mask_path = os.path.join(mask_dir, mask_name)
    if os.path.exists(mask_path):
        return mask_path
    return None

def calculate_cell_geometry(mask_path, image_path, output_dir, sample_name):
    os.makedirs(output_dir, exist_ok=True)
    mask = tifffile.imread(mask_path)
    image = tifffile.imread(image_path)
    if mask.shape != image.shape[1:]:
        raise ValueError(f"mask xy size does not match image: mask {mask.shape}, image {image.shape[1:]}")
    
    labels = np.unique(mask)
    labels = labels[labels != 0]
    
    cell_info = []
    for cell_label in labels:
        mask_single = (mask == cell_label)
        coords = np.argwhere(mask_single)
        y_min, x_min = coords.min(axis=0)
        cell_info.append((cell_label, y_min))
    
    cell_info.sort(key=lambda x: x[1])
    label_mapping = {}
    for i, (original_label, y_min) in enumerate(cell_info):
        label_mapping[original_label] = i + 1
    
    coordinate_mapping = {}
    csv_data = []
    
    for original_label, new_id in label_mapping.items():
        mask_single = (mask == original_label)
        coords = np.argwhere(mask_single)
        y_min, x_min = coords.min(axis=0)
        y_max, x_max = coords.max(axis=0)
        h = y_max - y_min + 1
        w = x_max - x_min + 1
        side = max(h, w)
        y_center = (y_min + y_max) // 2
        x_center = (x_min + x_max) // 2
        half = side // 2
        y1 = max(0, y_center - half)
        y2 = min(mask.shape[0], y1 + side)
        x1 = max(0, x_center - half)
        x2 = min(mask.shape[1], x1 + side)
        if y2 - y1 < side:
            y1 = max(0, y2 - side)
        if x2 - x1 < side:
            x1 = max(0, x2 - side)
        
        mask_crop = mask[y1:y2, x1:x2]
        image_crop = image[:, y1:y2, x1:x2]
        mask_single_crop = mask_crop == original_label
        
        # Estimate z-range from FITC signal
        z_projection = image_crop[:, mask_single_crop].sum(axis=1)  # sum of signal per z-slice
        threshold = 0.2 * np.max(z_projection)  # you can tune this
        z_indices = np.where(z_projection > threshold)[0]
        if len(z_indices) > 0:
            z_start, z_end = z_indices[0], z_indices[-1]
            cell_height_um = (z_end - z_start + 1) * z_step
            z_slices_count = z_end - z_start + 1
        else:
            cell_height_um = 0
            z_slices_count = 0
        
        # Extract 2D contour from mask and fit ellipse
        labeled_mask = label(mask_single_crop.astype(np.uint8))
        regions = regionprops(labeled_mask)
        if len(regions) > 0:
            region = regions[0]
            minor_axis = region.minor_axis_length * pixel_size_xy  # µm
            major_axis = region.major_axis_length * pixel_size_xy  # µm
            area_pixels = region.area
            area_um2 = area_pixels * (pixel_size_xy ** 2)
        else:
            minor_axis = 0
            major_axis = 0
            area_pixels = 0
            area_um2 = 0
        
        # Volume estimation - Triaxial ellipsoid (original method)
        volume_triaxial = 0
        if calculate_triaxial_ellipsoid:
            a = minor_axis / 2
            b = major_axis / 2
            c = cell_height_um / 2
            volume_triaxial = (4/3) * np.pi * a * b * c
        
        # Volume estimation - Prolate ellipsoid (rotation around major axis)
        volume_prolate = 0
        if calculate_prolate_ellipsoid:
            a = minor_axis / 2  # semi-minor axis (used twice)
            b = major_axis / 2  # semi-major axis (rotation axis)
            volume_prolate = (4/3) * np.pi * (a ** 2) * b
        
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
            },
            'geometry': {
                'height_um': float(cell_height_um),
                'major_axis_um': float(major_axis),
                'minor_axis_um': float(minor_axis),
                'volume_triaxial_ellipsoid_um3': float(volume_triaxial),
                'volume_prolate_ellipsoid_um3': float(volume_prolate)
            }
        }
        
        # Prepare CSV row
        csv_row = {
            'sample_name': sample_name,
            'cell_id': new_id,
            'original_mask_label': original_label,
            'height_um': round(cell_height_um, 4),
            'z_slices_count': z_slices_count,
            'major_axis_um': round(major_axis, 4),
            'minor_axis_um': round(minor_axis, 4),
            'area_um2': round(area_um2, 4),
            'triaxial_semi_a_um': round(minor_axis / 2, 4),
            'triaxial_semi_b_um': round(major_axis / 2, 4),
            'triaxial_semi_c_um': round(cell_height_um / 2, 4),
            'volume_triaxial_ellipsoid_um3': round(volume_triaxial, 6),
            'prolate_semi_a_um': round(minor_axis / 2, 4),
            'prolate_semi_b_um': round(major_axis / 2, 4),
            'volume_prolate_ellipsoid_um3': round(volume_prolate, 6),
            'pixel_size_xy': pixel_size_xy,
            'z_step_um': z_step,
            'z_threshold_factor': 0.2
        }
        csv_data.append(csv_row)
    
    # Save JSON
    mapping_file = os.path.join(output_dir, 'cell_geometry.json')
    with open(mapping_file, 'w') as f:
        json.dump(coordinate_mapping, f, indent=2)
    
    # Save CSV
    csv_file = os.path.join(output_dir, 'cell_volumes.csv')
    if csv_data:
        fieldnames = csv_data[0].keys()
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_data)
    
    return csv_data

if __name__ == '__main__':
    print(f"Pixel size XY: {pixel_size_xy} µm")
    print(f"Z-step size: {z_step} µm")
    print(f"Volume calculation modes:")
    print(f"  - Triaxial ellipsoid (z-height): {'Enabled' if calculate_triaxial_ellipsoid else 'Disabled'}")
    print(f"  - Prolate ellipsoid (major-axis rotation): {'Enabled' if calculate_prolate_ellipsoid else 'Disabled'}")
    print("Calculating cell geometry only (no image extraction)")
    
    all_csv_data = []
    
    for fname in os.listdir(deconv_dir):
        if not fname.lower().endswith('.tif') and not fname.lower().endswith('.tiff'):
            continue
        prefix, sidx = parse_deconv_filename(fname)
        if prefix is None or sidx is None:
            print(f"Ignoring: {fname}")
            continue
        mask_path = find_mask_file(prefix, sidx)
        if mask_path is None:
            print(f"No mask found: {prefix}, {sidx}")
            continue
        image_path = os.path.join(deconv_dir, fname)
        output_dir = os.path.join(output_base, f'{prefix}_{sidx}')
        sample_name = f'{prefix}_{sidx}'
        print(f"Processing: {fname} -> {os.path.basename(mask_path)} -> {output_dir}")
        try:
            csv_data = calculate_cell_geometry(mask_path, image_path, output_dir, sample_name)
            all_csv_data.extend(csv_data)
            print(f"Done: {output_dir} ({len(csv_data)} cells)")
        except Exception as e:
            print(f"Error processing {fname}: {e}")
    
    # Save combined CSV file
    if all_csv_data:
        combined_csv_file = os.path.join(output_base, 'all_cell_volumes.csv')
        fieldnames = all_csv_data[0].keys()
        with open(combined_csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_csv_data)
        print(f"\nCombined results saved to: {combined_csv_file}")
        print(f"Total cells processed: {len(all_csv_data)}") 