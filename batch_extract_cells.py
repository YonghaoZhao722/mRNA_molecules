import os
import tifffile
import numpy as np
import json
from scipy.ndimage import binary_dilation
from skimage.morphology import disk

# Configuration parameters
add_noise = False  # Set to False to disable Gaussian noise addition
base_dir = r'Y333 ATP6 ATP2'
deconv_dir = os.path.join(base_dir, 'binaries')
mask_dir = os.path.join(base_dir, 'aligned_masks')
output_base = os.path.join(base_dir, 'extracted_cells_binary_rm_noise')
os.makedirs(output_base, exist_ok=True)

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

def create_background_noise(image_slice, mask_single, band_width=5, mask=None):
    dilated_mask = binary_dilation(mask_single, structure=disk(band_width))
    band_region = dilated_mask & ~mask_single
    if np.sum(band_region) > 0:
        band_values = image_slice[band_region]
        mean_val = np.mean(band_values)
        std_val = np.std(band_values)
        noise_std = 0.5 * std_val
        return mean_val, noise_std
    else:
        if mask is not None:
            background_region = mask == 0
            if np.sum(background_region) > 0:
                bg_values = image_slice[background_region]
                mean_val = np.mean(bg_values)
                std_val = np.std(bg_values)
                return mean_val, 0.5 * std_val
        return 0, 0

def extract_cells(mask_path, image_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    mask = tifffile.imread(mask_path)
    image = tifffile.imread(image_path)
    if mask.shape != image.shape[1:]:
        raise ValueError(f"mask xy尺寸与原图不一致: mask {mask.shape}, image {image.shape[1:]}")
    
    labels = np.unique(mask)
    labels = labels[labels != 0]
    
    cell_info = []
    for label in labels:
        mask_single = (mask == label)
        coords = np.argwhere(mask_single)
        y_min, x_min = coords.min(axis=0)
        cell_info.append((label, y_min))
    
    # 按y_min坐标排序，然后分配从1开始的新编号
    cell_info.sort(key=lambda x: x[1])  # 按y_min排序
    label_mapping = {}
    for i, (original_label, y_min) in enumerate(cell_info):
        label_mapping[original_label] = i + 1  # 新编号从1开始
    
    coordinate_mapping = {}
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
        mask_crop = mask[y1:y2, x1:x2]
        image_crop = image[:, y1:y2, x1:x2]
        mask_single_crop = mask_crop == original_label
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
            
            result_slice[mask_single_crop] = current_slice[mask_single_crop]
            result[z] = result_slice
        output_path = os.path.join(output_dir, f'cell_{new_id:03d}.tif')
        tifffile.imwrite(output_path, result)
    mapping_file = os.path.join(output_dir, 'coordinate_mapping.json')
    with open(mapping_file, 'w') as f:
        json.dump(coordinate_mapping, f, indent=2)

if __name__ == '__main__':
    print(f"Gaussian noise addition: {'Enabled' if add_noise else 'Disabled'}")
    
    for fname in os.listdir(deconv_dir):
        if not fname.lower().endswith('.tif') and not fname.lower().endswith('.tiff'):
            continue
        prefix, sidx = parse_deconv_filename(fname)
        if prefix is None or sidx is None:
            print(f"跳过无法解析的文件: {fname}")
            continue
        mask_path = find_mask_file(prefix, sidx)
        if mask_path is None:
            print(f"未找到对应mask: {prefix}, {sidx}")
            continue
        image_path = os.path.join(deconv_dir, fname)
        output_dir = os.path.join(output_base, f'{prefix}_{sidx}')
        print(f"处理: {fname} -> {os.path.basename(mask_path)} 输出到 {output_dir}")
        try:
            extract_cells(mask_path, image_path, output_dir)
            print(f"完成: {output_dir}")
        except Exception as e:
            print(f"处理 {fname} 时出错: {e}")