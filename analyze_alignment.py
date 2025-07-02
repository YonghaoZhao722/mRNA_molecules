#!/usr/bin/env python3
"""
专门分析您的skeleton和spot文件对齐问题的脚本
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def load_skeleton_txt(file_path, mapping_file=None, cell_name="cell_001", pixel_size_xy=0.0645):
    """加载MitoGraph输出的skeleton文本文件，并转换为绝对坐标"""
    df = pd.read_csv(file_path, sep='\t')
    # 提取坐标列 (x, y, z) - 这些是相对于截取区域的坐标
    coords = df[['x', 'y', 'z']].values
    print(f"Skeleton数据加载完成: {len(coords)} 个点")
    print(f"相对坐标范围: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
    
    # 如果提供了mapping文件，转换为绝对坐标
    if mapping_file and os.path.exists(mapping_file):
        import json
        with open(mapping_file, 'r') as f:
            mapping = json.load(f)
        
        if cell_name in mapping:
            crop_info = mapping[cell_name]['crop_region']
            x_offset = crop_info['x_offset']
            y_offset = crop_info['y_offset']
            
            # 将像素偏移量转换为微米单位，并加到相对坐标上
            offset_x_um = x_offset * pixel_size_xy
            offset_y_um = y_offset * pixel_size_xy
            
            coords[:, 0] += offset_x_um  # X坐标加偏移
            coords[:, 1] += offset_y_um  # Y坐标加偏移
            # Z坐标不需要偏移，因为是3D图像的深度方向
            
            print(f"应用坐标偏移: X+{offset_x_um:.3f}μm (像素{x_offset}), Y+{offset_y_um:.3f}μm (像素{y_offset})")
            print(f"绝对坐标范围: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
        else:
            print(f"警告: 在mapping文件中未找到{cell_name}")
    else:
        if mapping_file:
            print(f"警告: mapping文件不存在: {mapping_file}")
        print("使用相对坐标（未应用偏移）")
    
    return coords

def load_spots_fishquant(file_path, cell_number=1, flip_y=True, mapping_data=None):
    """加载FISH-QUANT输出的spots文件，只提取指定细胞的数据"""
    # 读取文件找到SPOTS部分
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # 找到像素大小信息
    pixel_xy, pixel_z = None, None
    for i, line in enumerate(lines):
        if line.startswith('Pix-XY'):
            # 像素大小数值在下一行
            if i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                parts = next_line.split('\t')
                if len(parts) >= 2:
                    pixel_xy = float(parts[0])
                    pixel_z = float(parts[1])
            break
    
    print(f"FISH-QUANT像素大小: XY={pixel_xy}nm, Z={pixel_z}nm")
    
    # 找到指定细胞的SPOTS数据
    target_cell = f"CELL\tCell_{cell_number}"
    cell_found = False
    spots_start = -1
    spots_end = -1
    
    for i, line in enumerate(lines):
        if line.strip() == target_cell:
            cell_found = True
            print(f"找到目标细胞: Cell_{cell_number}")
            continue
        
        if cell_found and line.startswith('Pos_Y'):
            spots_start = i
            continue
            
        if cell_found and spots_start != -1 and line.startswith('CELL'):
            spots_end = i
            break
    
    # 如果没找到下一个CELL标记，说明是最后一个细胞
    if cell_found and spots_start != -1 and spots_end == -1:
        spots_end = len(lines)
    
    if not cell_found:
        raise ValueError(f"未找到Cell_{cell_number}的数据")
    
    if spots_start == -1:
        raise ValueError(f"未找到Cell_{cell_number}的SPOTS数据")
    
    # 读取该细胞的spots数据
    spots_data = []
    for i in range(spots_start + 1, spots_end):
        line = lines[i].strip()
        if line and not line.startswith('X_POS') and not line.startswith('Y_POS'):
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    y, x, z = float(parts[0]), float(parts[1]), float(parts[2])
                    spots_data.append([x, y, z])
                except ValueError:
                    continue
    
    coords = np.array(spots_data)
    
    # 如果需要翻转Y轴
    if flip_y and mapping_data:
        # 使用mapping数据中的crop区域信息进行精确翻转
        cell_name = f"cell_{cell_number:03d}"
        if cell_name in mapping_data:
            crop_info = mapping_data[cell_name]['crop_region']
            y_start = crop_info['y_start']  # 像素坐标
            y_end = crop_info['y_end']      # 像素坐标
            
            # 转换为纳米坐标 (像素 * 像素大小)
            y_start_nm = y_start * pixel_xy
            y_end_nm = y_end * pixel_xy
            
            # 基于crop区域进行翻转：new_y = (y_start + y_end) - old_y
            flip_center = y_start_nm + y_end_nm
            coords[:, 1] = flip_center - coords[:, 1]
            
            print(f"Y轴翻转: 基于crop区域 y_start={y_start_nm:.1f}nm, y_end={y_end_nm:.1f}nm")
        else:
            print(f"警告: 未找到{cell_name}的mapping信息，跳过Y轴翻转")
    elif flip_y:
        print("警告: 需要mapping数据进行精确Y轴翻转，使用简单翻转")
        # 简单翻转：基于当前细胞spots的范围
        if len(coords) > 0:
            y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
            flip_center = y_min + y_max
            coords[:, 1] = flip_center - coords[:, 1]
            print(f"Y轴翻转: 基于细胞范围，翻转中心={flip_center/2:.1f}nm")
    
    print(f"Cell_{cell_number} Spots数据加载完成: {len(coords)} 个点")
    print(f"坐标范围: X=[{coords[:,0].min():.1f}, {coords[:,0].max():.1f}], Y=[{coords[:,1].min():.1f}, {coords[:,1].max():.1f}], Z=[{coords[:,2].min():.1f}, {coords[:,2].max():.1f}]")
    
    return coords, pixel_xy, pixel_z

def analyze_unit_mismatch(skeleton_coords, spots_coords, pixel_xy, pixel_z):
    """分析单位不匹配问题"""
    print("\n=== 单位分析 ===")
    
    # 计算坐标范围
    skeleton_range = np.ptp(skeleton_coords, axis=0)
    spots_range = np.ptp(spots_coords, axis=0)
    
    print(f"坐标范围比较:")
    print(f"  Skeleton: X={skeleton_range[0]:.3f}μm, Y={skeleton_range[1]:.3f}μm, Z={skeleton_range[2]:.3f}μm")
    print(f"  Spots:    X={spots_range[0]:.1f}nm, Y={spots_range[1]:.1f}nm, Z={spots_range[2]:.1f}nm")
    
    # 转换spots到微米单位
    spots_coords_um = spots_coords / 1000.0
    spots_range_um = np.ptp(spots_coords_um, axis=0)
    
    print(f"  Spots(μm): X={spots_range_um[0]:.3f}μm, Y={spots_range_um[1]:.3f}μm, Z={spots_range_um[2]:.3f}μm")
    
    # 计算缩放比例
    scale_ratio = skeleton_range / spots_range_um
    print(f"  缩放比例: X={scale_ratio[0]:.3f}, Y={scale_ratio[1]:.3f}, Z={scale_ratio[2]:.3f}")
    
    # 检查pixel_size匹配
    expected_pixel_xy = pixel_xy / 1000.0  # 转换为微米
    expected_pixel_z = pixel_z / 1000.0
    
    print(f"\n像素大小检查:")
    print(f"  FISH-QUANT像素大小: XY={expected_pixel_xy:.4f}μm, Z={expected_pixel_z:.3f}μm")
    print(f"  您提供的像素大小: XY=0.0645μm, Z=0.2μm")
    
    if abs(expected_pixel_xy - 0.0645) < 0.001 and abs(expected_pixel_z - 0.2) < 0.01:
        print("  ✅ 像素大小匹配")
    else:
        print("  ⚠️ 像素大小不匹配")
    
    return spots_coords_um, scale_ratio

def calculate_alignment_after_correction(skeleton_coords, spots_coords_um):
    """计算校正后的对齐度量"""
    print("\n=== 对齐度量（微米单位）===")
    
    # 计算中心点
    skeleton_center = np.mean(skeleton_coords, axis=0)
    spots_center = np.mean(spots_coords_um, axis=0)
    
    print(f"中心点比较:")
    print(f"  Skeleton: X={skeleton_center[0]:.3f}, Y={skeleton_center[1]:.3f}, Z={skeleton_center[2]:.3f}")
    print(f"  Spots:    X={spots_center[0]:.3f}, Y={spots_center[1]:.3f}, Z={spots_center[2]:.3f}")
    
    center_diff = skeleton_center - spots_center
    center_distance = np.linalg.norm(center_diff)
    print(f"  中心距离: {center_distance:.3f}μm")
    
    # 计算最近邻距离（简化版本，使用采样）
    if len(skeleton_coords) > 1000:
        skeleton_sample = skeleton_coords[::len(skeleton_coords)//1000]
    else:
        skeleton_sample = skeleton_coords
    
    if len(spots_coords_um) > 500:
        spots_sample = spots_coords_um[::len(spots_coords_um)//500]
    else:
        spots_sample = spots_coords_um
    
    from scipy.spatial.distance import cdist
    distances = cdist(spots_sample, skeleton_sample)
    min_distances = np.min(distances, axis=1)
    
    print(f"\n最近邻距离统计:")
    print(f"  平均距离: {np.mean(min_distances):.3f}μm")
    print(f"  中位数距离: {np.median(min_distances):.3f}μm")
    print(f"  最大距离: {np.max(min_distances):.3f}μm")
    print(f"  1μm内的spots: {np.sum(min_distances < 1.0) / len(min_distances) * 100:.1f}%")
    print(f"  0.5μm内的spots: {np.sum(min_distances < 0.5) / len(min_distances) * 100:.1f}%")
    print(f"  0.1μm内的spots: {np.sum(min_distances < 0.1) / len(min_distances) * 100:.1f}%")
    
    return center_diff, min_distances

def plot_alignment_analysis(skeleton_coords, spots_coords_um, min_distances):
    """绘制对齐分析图"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    
    # XY投影 - 设置等比例
    axes[0, 0].scatter(skeleton_coords[:, 0], skeleton_coords[:, 1], 
                      c='red', alpha=0.6, s=1, label='Skeleton')
    axes[0, 0].scatter(spots_coords_um[:, 0], spots_coords_um[:, 1], 
                      c='blue', alpha=0.8, s=20, label='Spots')
    axes[0, 0].set_xlabel('X (μm)')
    axes[0, 0].set_ylabel('Y (μm)')
    axes[0, 0].set_title('XY Projection')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_aspect('equal', adjustable='box')  # 设置1:1比例
    
    # XZ投影 - 设置等比例
    axes[0, 1].scatter(skeleton_coords[:, 0], skeleton_coords[:, 2], 
                      c='red', alpha=0.6, s=1, label='Skeleton')
    axes[0, 1].scatter(spots_coords_um[:, 0], spots_coords_um[:, 2], 
                      c='blue', alpha=0.8, s=20, label='Spots')
    axes[0, 1].set_xlabel('X (μm)')
    axes[0, 1].set_ylabel('Z (μm)')
    axes[0, 1].set_title('XZ Projection')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_aspect('equal', adjustable='box')  # 设置1:1比例
    
    # YZ投影 - 设置等比例
    axes[1, 0].scatter(skeleton_coords[:, 1], skeleton_coords[:, 2], 
                      c='red', alpha=0.6, s=1, label='Skeleton')
    axes[1, 0].scatter(spots_coords_um[:, 1], spots_coords_um[:, 2], 
                      c='blue', alpha=0.8, s=20, label='Spots')
    axes[1, 0].set_xlabel('Y (μm)')
    axes[1, 0].set_ylabel('Z (μm)')
    axes[1, 0].set_title('YZ Projection')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_aspect('equal', adjustable='box')  # 设置1:1比例
    
    # 距离分布直方图
    axes[1, 1].hist(min_distances, bins=30, alpha=0.7, edgecolor='black')
    axes[1, 1].set_xlabel('Distance to Nearest Skeleton Point (μm)')
    axes[1, 1].set_ylabel('Number of Spots')
    axes[1, 1].set_title('Distance Distribution')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axvline(np.mean(min_distances), color='red', linestyle='--', 
                      label=f'Mean: {np.mean(min_distances):.3f}μm')
    axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig('alignment_analysis.png', dpi=150, bbox_inches='tight')
    print(f"\n图表已保存到: alignment_analysis.png")

def main():
    import sys
    
    # 支持命令行参数指定细胞编号
    if len(sys.argv) >= 2:
        cell_num = sys.argv[1].zfill(3)  # 补齐为3位数字，如 3 -> 003
        skeleton_file = f"Y333 ATP6 ATP2/extracted_cells/yWL333_cy3_ATP2_cy5_ATP6MS2_1_s1/cell_{cell_num}.txt"
    else:
        # 默认分析cell_003
        skeleton_file = "Y333 ATP6 ATP2/extracted_cells/yWL333_cy3_ATP2_cy5_ATP6MS2_1_s1/cell_003.txt"
    
    spots_file = "Y333 ATP6 ATP2/atp6_spots/yWL333_cy3_ATP2_cy5_ATP6MS2_1_DIC_s1_spots_250622.txt"
    
    print(f"分析细胞: {os.path.basename(skeleton_file)}")
    print(f"使用方法: python analyze_alignment.py [细胞编号]")
    print(f"例如: python analyze_alignment.py 1  # 分析cell_001")
    print(f"注意: 默认启用Y轴翻转以修正坐标系统差异")
    
    print("=== Skeleton和Spot坐标对齐分析 ===")
    print(f"Skeleton文件: {skeleton_file}")
    print(f"Spots文件: {spots_file}")
    
    try:
        # 构建mapping文件路径
        skeleton_dir = os.path.dirname(skeleton_file)
        mapping_file = os.path.join(skeleton_dir, 'coordinate_mapping.json')
        
        # 加载mapping数据
        mapping_data = None
        if os.path.exists(mapping_file):
            import json
            with open(mapping_file, 'r') as f:
                mapping_data = json.load(f)
        
        # 从skeleton文件名自动识别细胞编号
        skeleton_basename = os.path.basename(skeleton_file)
        if skeleton_basename.startswith('cell_'):
            # 提取细胞编号，如 cell_001.txt -> 1, cell_003.txt -> 3
            cell_id_str = skeleton_basename.split('_')[1].split('.')[0]
            cell_number = int(cell_id_str.lstrip('0')) if cell_id_str != '000' else 0
            cell_name = f"cell_{cell_id_str}"
        else:
            cell_number = 1
            cell_name = "cell_001"
        
        print(f"自动识别细胞: {cell_name} (Cell_{cell_number})")
        
        # 加载数据 - 先加载spots获取像素大小，使用mapping数据进行精确Y轴翻转
        spots_coords, pixel_xy, pixel_z = load_spots_fishquant(spots_file, cell_number=cell_number, flip_y=True, mapping_data=mapping_data)
        # 然后加载skeleton并应用坐标转换
        skeleton_coords = load_skeleton_txt(skeleton_file, mapping_file=mapping_file, 
                                          cell_name=cell_name, pixel_size_xy=pixel_xy/1000)
        
        # 分析单位不匹配
        spots_coords_um, scale_ratio = analyze_unit_mismatch(skeleton_coords, spots_coords, pixel_xy, pixel_z)
        
        # 计算对齐度量
        center_diff, min_distances = calculate_alignment_after_correction(skeleton_coords, spots_coords_um)
        
        # 绘制分析图
        plot_alignment_analysis(skeleton_coords, spots_coords_um, min_distances)
        
        # 结论和建议
        print(f"\n=== 分析结论 ===")
        print(f"1. 单位问题: FISH-QUANT输出的是纳米单位，MitoGraph输出的是微米单位")
        print(f"2. 转换方法: spots坐标除以1000即可转换为微米单位")
        print(f"3. 转换后中心距离: {np.linalg.norm(center_diff):.3f}μm")
        
        if np.linalg.norm(center_diff) < 0.5:
            print(f"4. ✅ 转换后坐标对齐良好")
        else:
            print(f"4. ⚠️ 转换后仍有偏移，可能需要进一步校正")
            print(f"   建议偏移量: X={-center_diff[0]:.3f}μm, Y={-center_diff[1]:.3f}μm, Z={-center_diff[2]:.3f}μm")
        
        # 生成校正后的spots文件
        print(f"\n=== 生成校正文件 ===")
        
        # 创建校正后的spots坐标（微米单位）
        spots_corrected = spots_coords_um - center_diff
        
        # 保存为简单的CSV文件
        corrected_df = pd.DataFrame(spots_corrected, columns=['x_um', 'y_um', 'z_um'])
        corrected_df.to_csv('spots_corrected.csv', index=False)
        print(f"校正后的spots坐标已保存到: spots_corrected.csv")
        
        # 显示最终对齐统计
        from scipy.spatial.distance import cdist
        final_distances = cdist(spots_corrected[:500], skeleton_coords[:1000])
        final_min_distances = np.min(final_distances, axis=1)
        
        print(f"\n=== 最终对齐统计 ===")
        print(f"平均距离: {np.mean(final_min_distances):.3f}μm")
        print(f"1μm内的spots: {np.sum(final_min_distances < 1.0) / len(final_min_distances) * 100:.1f}%")
        print(f"0.5μm内的spots: {np.sum(final_min_distances < 0.5) / len(final_min_distances) * 100:.1f}%")
        print(f"0.1μm内的spots: {np.sum(final_min_distances < 0.1) / len(final_min_distances) * 100:.1f}%")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 