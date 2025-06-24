import os
import tifffile
import numpy as np
import json
from scipy.ndimage import binary_dilation, binary_erosion
from skimage.morphology import disk

# 路径
mask_path = r'sample_atp3/committed_objects_shifted.tif'
image_path = r'sample_atp3/atp3RL.tif'
output_dir = r'sample_atp3/mitograph'
os.makedirs(output_dir, exist_ok=True)

# 读取mask和原图
mask = tifffile.imread(mask_path)
image = tifffile.imread(image_path)

# 检查mask和image的xy尺寸
if mask.shape != image.shape[1:]:
    raise ValueError(f"mask xy尺寸与原图不一致: mask {mask.shape}, image {image.shape[1:]}")

labels = np.unique(mask)
labels = labels[labels != 0]  # 去除背景

print(f"检测到{len(labels)}个细胞")
print(f"细胞标签: {labels}")

# 保存坐标映射信息
coordinate_mapping = {}

def create_background_noise(image_slice, mask_single, band_width=5):
    """
    创建背景噪声，模拟ImageJ宏中的逻辑
    在mask周围创建band区域，计算统计信息并生成噪声
    """
    # 创建扩展的mask (dilated)
    dilated_mask = binary_dilation(mask_single, structure=disk(band_width))
    # band区域 = dilated - original
    band_region = dilated_mask & ~mask_single
    
    if np.sum(band_region) > 0:
        # 计算band区域的统计信息
        band_values = image_slice[band_region]
        mean_val = np.mean(band_values)
        std_val = np.std(band_values)
        
        # 生成噪声 (类似ImageJ中的0.5*std)
        noise_std = 0.5 * std_val
        return mean_val, noise_std
    else:
        # 如果没有band区域，使用整个图像的背景统计
        background_region = mask == 0
        if np.sum(background_region) > 0:
            bg_values = image_slice[background_region]
            mean_val = np.mean(bg_values)
            std_val = np.std(bg_values)
            return mean_val, 0.5 * std_val
        else:
            return 0, 0

for label in labels:
    print(f"处理细胞 {label}...")
    
    mask_single = (mask == label)
    coords = np.argwhere(mask_single)
    y_min, x_min = coords.min(axis=0)
    y_max, x_max = coords.max(axis=0)
    
    # 计算正方形边长
    h = y_max - y_min + 1
    w = x_max - x_min + 1
    side = max(h, w)
    
    # 以中心为基准扩展为正方形
    y_center = (y_min + y_max) // 2
    x_center = (x_min + x_max) // 2
    half = side // 2
    
    y1 = max(0, y_center - half)
    y2 = min(mask.shape[0], y1 + side)
    x1 = max(0, x_center - half)
    x2 = min(mask.shape[1], x1 + side)
    
    # 修正边界（防止超出图像）
    if y2 - y1 < side:
        y1 = max(0, y2 - side)
    if x2 - x1 < side:
        x1 = max(0, x2 - side)
    
    # 保存坐标映射信息
    coordinate_mapping[f'cell_{label:03d}'] = {
        'cell_id': int(label),
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
    
    # 裁剪
    mask_crop = mask[y1:y2, x1:x2]
    image_crop = image[:, y1:y2, x1:x2]
    mask_single_crop = mask_crop == label
    
    # 创建结果数组
    result = np.zeros_like(image_crop, dtype=np.uint16)
    
    # 对每个z层处理
    for z in range(image.shape[0]):
        current_slice = image_crop[z]
        result_slice = np.zeros_like(current_slice, dtype=np.uint16)
        
        # 计算背景噪声参数
        mean_noise, std_noise = create_background_noise(current_slice, mask_single_crop)
        
        # 生成噪声填充
        noise_shape = current_slice.shape
        noise = np.random.normal(mean_noise, std_noise, noise_shape)
        noise = np.clip(noise, 0, 65535).astype(np.uint16)  # 限制在uint16范围内
        
        # 填充结果：
        # 1. 首先用噪声填充整个区域
        result_slice[:] = noise
        # 2. 然后在mask区域覆盖真实数据
        result_slice[mask_single_crop] = current_slice[mask_single_crop]
        
        result[z] = result_slice
    
    output_path = os.path.join(output_dir, f'cell_{label:03d}.tif')
    tifffile.imwrite(output_path, result)
    print(f"已保存: {output_path}")

# 保存坐标映射文件
mapping_file = os.path.join(output_dir, 'coordinate_mapping.json')
with open(mapping_file, 'w') as f:
    json.dump(coordinate_mapping, f, indent=2)
print(f"坐标映射文件已保存: {mapping_file}")

print("全部完成！")