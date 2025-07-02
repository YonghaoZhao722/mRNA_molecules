#!/usr/bin/env python3
"""
Coordinate Alignment Diagnostic Tool
用于诊断MitoGraph skeleton和FISH-QUANT spot之间的坐标对齐问题

Usage:
python coordinate_alignment_diagnostic.py --skeleton skeleton.vtk --spots spots.csv --pixel-xy 0.0645 --pixel-z 0.2
"""

import argparse
import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist


class CoordinateAlignmentDiagnostic:
    def __init__(self, pixel_size_xy=0.0645, pixel_size_z=0.2):
        """
        初始化坐标对齐诊断工具
        
        Args:
            pixel_size_xy: XY方向像素大小 (微米)
            pixel_size_z: Z方向像素大小 (微米)
        """
        self.pixel_size_xy = pixel_size_xy
        self.pixel_size_z = pixel_size_z
        
    def load_skeleton_vtk(self, vtk_file):
        """加载VTK格式的skeleton数据"""
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(vtk_file)
        reader.Update()
        
        polydata = reader.GetOutput()
        points = polydata.GetPoints()
        
        # 转换为numpy数组
        skeleton_coords = vtk_to_numpy(points.GetData())
        
        print(f"Skeleton坐标范围:")
        print(f"  X: {skeleton_coords[:, 0].min():.3f} - {skeleton_coords[:, 0].max():.3f}")
        print(f"  Y: {skeleton_coords[:, 1].min():.3f} - {skeleton_coords[:, 1].max():.3f}")
        print(f"  Z: {skeleton_coords[:, 2].min():.3f} - {skeleton_coords[:, 2].max():.3f}")
        
        return skeleton_coords
    
    def load_spots_csv(self, csv_file):
        """加载CSV格式的spot数据"""
        # 尝试不同的分隔符
        for sep in [',', '\t', ';']:
            try:
                spots_df = pd.read_csv(csv_file, sep=sep)
                if len(spots_df.columns) >= 3:
                    break
            except:
                continue
        
        # 查找坐标列
        coord_cols = []
        for possible_names in [['x', 'y', 'z'], ['X', 'Y', 'Z'], ['pos_x', 'pos_y', 'pos_z']]:
            if all(col in spots_df.columns for col in possible_names):
                coord_cols = possible_names
                break
        
        if not coord_cols:
            # 假设前三列是坐标
            coord_cols = spots_df.columns[:3].tolist()
            print(f"警告: 自动选择前三列作为坐标: {coord_cols}")
        
        spots_coords = spots_df[coord_cols].values
        
        print(f"Spot坐标范围:")
        print(f"  X: {spots_coords[:, 0].min():.3f} - {spots_coords[:, 0].max():.3f}")
        print(f"  Y: {spots_coords[:, 1].min():.3f} - {spots_coords[:, 1].max():.3f}")
        print(f"  Z: {spots_coords[:, 2].min():.3f} - {spots_coords[:, 2].max():.3f}")
        
        return spots_coords, coord_cols
    
    def analyze_coordinate_scaling(self, skeleton_coords, spots_coords):
        """分析坐标缩放差异"""
        print("\n=== 坐标缩放分析 ===")
        
        # 计算坐标范围比例
        skeleton_range = np.ptp(skeleton_coords, axis=0)
        spots_range = np.ptp(spots_coords, axis=0)
        
        scale_ratio = skeleton_range / spots_range
        
        print(f"坐标范围:")
        print(f"  Skeleton: X={skeleton_range[0]:.3f}, Y={skeleton_range[1]:.3f}, Z={skeleton_range[2]:.3f}")
        print(f"  Spots:    X={spots_range[0]:.3f}, Y={spots_range[1]:.3f}, Z={spots_range[2]:.3f}")
        print(f"  比例:     X={scale_ratio[0]:.3f}, Y={scale_ratio[1]:.3f}, Z={scale_ratio[2]:.3f}")
        
        # 检查是否存在单位转换问题
        expected_ratio_nm_to_um = [1000, 1000, 1000]  # 如果一个是nm，另一个是μm
        expected_ratio_pixel = [self.pixel_size_xy, self.pixel_size_xy, self.pixel_size_z]
        
        print(f"\n可能的单位问题:")
        print(f"  如果存在nm/μm转换: {np.allclose(scale_ratio, expected_ratio_nm_to_um, rtol=0.1)}")
        print(f"  如果存在像素/物理单位: {np.allclose(scale_ratio, expected_ratio_pixel, rtol=0.1)}")
        print(f"  如果存在像素/物理单位(倒数): {np.allclose(scale_ratio, 1/np.array(expected_ratio_pixel), rtol=0.1)}")
        
        return scale_ratio
    
    def suggest_corrections(self, skeleton_coords, spots_coords, scale_ratio):
        """建议校正方案"""
        print("\n=== 校正建议 ===")
        
        # 方案1: 缩放skeleton匹配spots
        skeleton_scaled_1 = skeleton_coords / scale_ratio
        
        # 方案2: 缩放spots匹配skeleton  
        spots_scaled_2 = spots_coords * scale_ratio
        
        # 方案3: 单位转换
        if np.allclose(scale_ratio, [1000, 1000, 1000], rtol=0.1):
            print("方案A: Skeleton可能是nm单位，spots是μm单位")
            skeleton_corrected = skeleton_coords / 1000
        elif np.allclose(scale_ratio, [0.001, 0.001, 0.001], rtol=0.1):
            print("方案B: Spots可能是nm单位，skeleton是μm单位")
            spots_corrected = spots_coords / 1000
        
        # 方案4: 像素单位转换
        if np.allclose(scale_ratio[0], 1/self.pixel_size_xy, rtol=0.1):
            print("方案C: Spots可能是像素单位，需要转换为物理单位")
            spots_corrected = spots_coords * [self.pixel_size_xy, self.pixel_size_xy, self.pixel_size_z]
        elif np.allclose(scale_ratio[0], self.pixel_size_xy, rtol=0.1):
            print("方案D: Skeleton可能需要从物理单位转换回像素单位")
            skeleton_corrected = skeleton_coords / [self.pixel_size_xy, self.pixel_size_xy, self.pixel_size_z]
        
        return {
            'skeleton_scaled_to_spots': skeleton_scaled_1,
            'spots_scaled_to_skeleton': spots_scaled_2
        }
    
    def calculate_alignment_metrics(self, skeleton_coords, spots_coords):
        """计算对齐度量"""
        print("\n=== 对齐度量 ===")
        
        # 计算最近邻距离
        distances = cdist(spots_coords, skeleton_coords)
        min_distances = np.min(distances, axis=1)
        
        metrics = {
            'mean_distance': np.mean(min_distances),
            'median_distance': np.median(min_distances),
            'std_distance': np.std(min_distances),
            'max_distance': np.max(min_distances),
            'spots_within_1um': np.sum(min_distances < 1.0) / len(min_distances) * 100,
            'spots_within_0.5um': np.sum(min_distances < 0.5) / len(min_distances) * 100,
            'spots_within_0.1um': np.sum(min_distances < 0.1) / len(min_distances) * 100
        }
        
        print(f"平均距离: {metrics['mean_distance']:.3f} μm")
        print(f"中位数距离: {metrics['median_distance']:.3f} μm")
        print(f"最大距离: {metrics['max_distance']:.3f} μm")
        print(f"1μm内的spots: {metrics['spots_within_1um']:.1f}%")
        print(f"0.5μm内的spots: {metrics['spots_within_0.5um']:.1f}%")
        print(f"0.1μm内的spots: {metrics['spots_within_0.1um']:.1f}%")
        
        return metrics, min_distances
    
    def plot_alignment(self, skeleton_coords, spots_coords, min_distances, output_file=None):
        """绘制对齐分析图"""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # XY投影
        axes[0, 0].scatter(skeleton_coords[:, 0], skeleton_coords[:, 1], 
                          c='red', alpha=0.6, s=1, label='Skeleton')
        axes[0, 0].scatter(spots_coords[:, 0], spots_coords[:, 1], 
                          c='blue', alpha=0.8, s=20, label='Spots')
        axes[0, 0].set_xlabel('X (μm)')
        axes[0, 0].set_ylabel('Y (μm)')
        axes[0, 0].set_title('XY Projection')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # XZ投影
        axes[0, 1].scatter(skeleton_coords[:, 0], skeleton_coords[:, 2], 
                          c='red', alpha=0.6, s=1, label='Skeleton')
        axes[0, 1].scatter(spots_coords[:, 0], spots_coords[:, 2], 
                          c='blue', alpha=0.8, s=20, label='Spots')
        axes[0, 1].set_xlabel('X (μm)')
        axes[0, 1].set_ylabel('Z (μm)')
        axes[0, 1].set_title('XZ Projection')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 距离分布直方图
        axes[1, 0].hist(min_distances, bins=50, alpha=0.7, edgecolor='black')
        axes[1, 0].set_xlabel('Distance to Nearest Skeleton Point (μm)')
        axes[1, 0].set_ylabel('Number of Spots')
        axes[1, 0].set_title('Distance Distribution')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].axvline(np.mean(min_distances), color='red', linestyle='--', 
                          label=f'Mean: {np.mean(min_distances):.3f}')
        axes[1, 0].legend()
        
        # 3D散点图(降维到2D显示)
        # 使用颜色编码Z坐标
        scatter = axes[1, 1].scatter(spots_coords[:, 0], spots_coords[:, 1], 
                                    c=min_distances, cmap='viridis', s=20, alpha=0.8)
        axes[1, 1].scatter(skeleton_coords[:, 0], skeleton_coords[:, 1], 
                          c='red', alpha=0.3, s=1)
        axes[1, 1].set_xlabel('X (μm)')
        axes[1, 1].set_ylabel('Y (μm)')
        axes[1, 1].set_title('Spots Colored by Distance to Skeleton')
        plt.colorbar(scatter, ax=axes[1, 1], label='Distance (μm)')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"图表已保存到: {output_file}")
        else:
            plt.show()
    
    def generate_correction_script(self, correction_type, scale_factor=None):
        """生成校正脚本"""
        script_content = f"""#!/usr/bin/env python3
# Coordinate Correction Script
# Generated by CoordinateAlignmentDiagnostic

import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

def correct_skeleton_coordinates(vtk_file, output_file, correction_type="{correction_type}"):
    # 读取skeleton
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_file)
    reader.Update()
    
    polydata = reader.GetOutput()
    points = polydata.GetPoints()
    coords = vtk_to_numpy(points.GetData())
    
    # 应用校正
    if correction_type == "scale":
        scale_factor = {scale_factor if scale_factor else [1.0, 1.0, 1.0]}
        corrected_coords = coords / scale_factor
    elif correction_type == "nm_to_um":
        corrected_coords = coords / 1000.0
    elif correction_type == "pixel_to_physical":
        pixel_size = [{self.pixel_size_xy}, {self.pixel_size_xy}, {self.pixel_size_z}]
        corrected_coords = coords * pixel_size
    
    # 更新坐标
    corrected_points = numpy_to_vtk(corrected_coords)
    new_points = vtk.vtkPoints()
    new_points.SetData(corrected_points)
    polydata.SetPoints(new_points)
    
    # 保存校正后的文件
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(polydata)
    writer.Write()
    
    print(f"校正后的skeleton已保存到: {{output_file}}")

def correct_spots_coordinates(csv_file, output_file, correction_type="{correction_type}"):
    # 读取spots
    spots_df = pd.read_csv(csv_file)
    
    # 应用校正 (这里需要根据具体的CSV格式调整)
    # ...
    
    spots_df.to_csv(output_file, index=False)
    print(f"校正后的spots已保存到: {{output_file}}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python correction_script.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if input_file.endswith('.vtk'):
        correct_skeleton_coordinates(input_file, output_file)
    else:
        correct_spots_coordinates(input_file, output_file)
"""
        
        with open('coordinate_correction_script.py', 'w') as f:
            f.write(script_content)
        
        print("校正脚本已生成: coordinate_correction_script.py")


def main():
    parser = argparse.ArgumentParser(description='诊断skeleton和spot坐标对齐问题')
    parser.add_argument('--skeleton', required=True, help='Skeleton VTK文件路径')
    parser.add_argument('--spots', required=True, help='Spots CSV文件路径')
    parser.add_argument('--pixel-xy', type=float, default=0.0645, help='XY像素大小 (μm)')
    parser.add_argument('--pixel-z', type=float, default=0.2, help='Z像素大小 (μm)')
    parser.add_argument('--output-plot', help='输出图表文件路径')
    parser.add_argument('--generate-correction', action='store_true', help='生成校正脚本')
    
    args = parser.parse_args()
    
    # 初始化诊断工具
    diagnostic = CoordinateAlignmentDiagnostic(args.pixel_xy, args.pixel_z)
    
    print("=== 坐标对齐诊断工具 ===")
    print(f"像素大小: XY={args.pixel_xy}μm, Z={args.pixel_z}μm")
    
    # 加载数据
    print(f"\n加载skeleton: {args.skeleton}")
    skeleton_coords = diagnostic.load_skeleton_vtk(args.skeleton)
    
    print(f"\n加载spots: {args.spots}")
    spots_coords, coord_cols = diagnostic.load_spots_csv(args.spots)
    
    # 分析坐标缩放
    scale_ratio = diagnostic.analyze_coordinate_scaling(skeleton_coords, spots_coords)
    
    # 计算对齐度量
    metrics, min_distances = diagnostic.calculate_alignment_metrics(skeleton_coords, spots_coords)
    
    # 建议校正
    corrections = diagnostic.suggest_corrections(skeleton_coords, spots_coords, scale_ratio)
    
    # 绘制分析图
    diagnostic.plot_alignment(skeleton_coords, spots_coords, min_distances, args.output_plot)
    
    # 生成校正脚本
    if args.generate_correction:
        diagnostic.generate_correction_script("scale", scale_ratio)


if __name__ == "__main__":
    main() 