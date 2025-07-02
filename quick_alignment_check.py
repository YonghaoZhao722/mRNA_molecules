#!/usr/bin/env python3
"""
Quick Alignment Check
快速检查skeleton和spot的对齐情况

使用方法:
python quick_alignment_check.py /path/to/skeleton.vtk /path/to/spots.csv
"""

import sys
import os
import numpy as np
import pandas as pd

def load_vtk_simple(vtk_file):
    """简化的VTK文件读取（不依赖VTK库）"""
    coords = []
    with open(vtk_file, 'r') as f:
        in_points = False
        point_count = 0
        read_count = 0
        
        for line in f:
            line = line.strip()
            
            if line.startswith('POINTS'):
                point_count = int(line.split()[1])
                in_points = True
                continue
            
            if in_points and read_count < point_count:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                        coords.append([x, y, z])
                        read_count += 1
                    except ValueError:
                        continue
            
            if read_count >= point_count:
                break
    
    return np.array(coords)

def load_csv_simple(csv_file):
    """简化的CSV文件读取"""
    # 尝试不同的分隔符
    for sep in [',', '\t', ';']:
        try:
            df = pd.read_csv(csv_file, sep=sep)
            if len(df.columns) >= 3:
                break
        except:
            continue
    
    # 尝试找到坐标列
    coord_cols = []
    for possible in [['x', 'y', 'z'], ['X', 'Y', 'Z'], ['pos_x', 'pos_y', 'pos_z']]:
        if all(col in df.columns for col in possible):
            coord_cols = possible
            break
    
    if not coord_cols:
        # 使用前三列
        coord_cols = df.columns[:3].tolist()
        print(f"自动选择坐标列: {coord_cols}")
    
    coords = df[coord_cols].values
    return coords, coord_cols

def quick_alignment_analysis(skeleton_coords, spots_coords):
    """快速对齐分析"""
    print("\n=== 快速对齐分析 ===")
    
    # 坐标范围
    skeleton_range = np.ptp(skeleton_coords, axis=0)
    spots_range = np.ptp(spots_coords, axis=0)
    
    print(f"坐标范围比较:")
    print(f"  Skeleton - X: {skeleton_range[0]:.3f}, Y: {skeleton_range[1]:.3f}, Z: {skeleton_range[2]:.3f}")
    print(f"  Spots    - X: {spots_range[0]:.3f}, Y: {spots_range[1]:.3f}, Z: {spots_range[2]:.3f}")
    
    # 缩放比例
    scale_ratio = skeleton_range / spots_range
    print(f"  缩放比例 - X: {scale_ratio[0]:.3f}, Y: {scale_ratio[1]:.3f}, Z: {scale_ratio[2]:.3f}")
    
    # 中心点比较
    skeleton_center = np.mean(skeleton_coords, axis=0)
    spots_center = np.mean(spots_coords, axis=0)
    
    print(f"\n中心点比较:")
    print(f"  Skeleton - X: {skeleton_center[0]:.3f}, Y: {skeleton_center[1]:.3f}, Z: {skeleton_center[2]:.3f}")
    print(f"  Spots    - X: {spots_center[0]:.3f}, Y: {spots_center[1]:.3f}, Z: {spots_center[2]:.3f}")
    
    center_diff = skeleton_center - spots_center
    print(f"  中心偏移 - X: {center_diff[0]:.3f}, Y: {center_diff[1]:.3f}, Z: {center_diff[2]:.3f}")
    
    # 可能的问题诊断
    print(f"\n=== 问题诊断 ===")
    
    # 检查单位问题
    if np.allclose(scale_ratio, [1000, 1000, 1000], rtol=0.1):
        print("⚠️ 可能的问题: Skeleton是nm单位，spots是μm单位")
        print("   建议: 将skeleton坐标除以1000")
    elif np.allclose(scale_ratio, [0.001, 0.001, 0.001], rtol=0.1):
        print("⚠️ 可能的问题: Spots是nm单位，skeleton是μm单位") 
        print("   建议: 将spots坐标除以1000")
    
    # 检查像素单位问题 (基于您的像素大小)
    pixel_xy, pixel_z = 0.0645, 0.2
    if np.allclose(scale_ratio[0], 1/pixel_xy, rtol=0.2):
        print("⚠️ 可能的问题: Spots是像素单位，skeleton是物理单位")
        print(f"   建议: 将spots坐标乘以 [{pixel_xy}, {pixel_xy}, {pixel_z}]")
    elif np.allclose(scale_ratio[0], pixel_xy, rtol=0.2):
        print("⚠️ 可能的问题: Skeleton应该转换为像素单位")
        print(f"   建议: 将skeleton坐标除以 [{pixel_xy}, {pixel_xy}, {pixel_z}]")
    
    # 检查严重的尺寸差异
    if any(abs(scale_ratio - 1) > 0.5):
        print("⚠️ 检测到显著的尺寸差异，可能存在单位或缩放问题")
    
    # 检查中心偏移
    if np.linalg.norm(center_diff) > 1.0:  # 如果中心偏移超过1μm
        print("⚠️ 检测到显著的中心偏移，可能存在坐标系不匹配")
    
    return scale_ratio, center_diff

def suggest_fix(scale_ratio):
    """建议修复方案"""
    print(f"\n=== 修复建议 ===")
    
    # 基于缩放比例建议修复
    if np.allclose(scale_ratio, [1, 1, 1], rtol=0.1):
        print("✅ 坐标尺寸匹配良好，可能只需要微调对齐")
    elif np.allclose(scale_ratio, [1000, 1000, 1000], rtol=0.1):
        print("🔧 建议修复: skeleton_corrected = skeleton_coords / 1000")
    elif np.allclose(scale_ratio, [0.001, 0.001, 0.001], rtol=0.1):
        print("🔧 建议修复: spots_corrected = spots_coords / 1000")
    else:
        print(f"🔧 建议修复: skeleton_corrected = skeleton_coords / {scale_ratio}")
        print(f"   或者: spots_corrected = spots_coords * {scale_ratio}")

def main():
    if len(sys.argv) != 3:
        print("使用方法: python quick_alignment_check.py skeleton.vtk spots.csv")
        sys.exit(1)
    
    skeleton_file = sys.argv[1]
    spots_file = sys.argv[2]
    
    # 检查文件存在
    if not os.path.exists(skeleton_file):
        print(f"错误: 找不到skeleton文件 {skeleton_file}")
        sys.exit(1)
    
    if not os.path.exists(spots_file):
        print(f"错误: 找不到spots文件 {spots_file}")
        sys.exit(1)
    
    print("=== 快速坐标对齐检查 ===")
    print(f"Skeleton文件: {skeleton_file}")
    print(f"Spots文件: {spots_file}")
    
    try:
        # 加载数据
        print(f"\n正在加载skeleton数据...")
        skeleton_coords = load_vtk_simple(skeleton_file)
        print(f"加载了 {len(skeleton_coords)} 个skeleton点")
        
        print(f"\n正在加载spots数据...")
        spots_coords, coord_cols = load_csv_simple(spots_file)
        print(f"加载了 {len(spots_coords)} 个spots，坐标列: {coord_cols}")
        
        # 分析对齐
        scale_ratio, center_diff = quick_alignment_analysis(skeleton_coords, spots_coords)
        
        # 建议修复
        suggest_fix(scale_ratio)
        
        print(f"\n✅ 检查完成")
        print(f"如需详细分析，请使用: python coordinate_alignment_diagnostic.py --skeleton {skeleton_file} --spots {spots_file}")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 