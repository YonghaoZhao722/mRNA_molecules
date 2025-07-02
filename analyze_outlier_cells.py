import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json

def check_cells_with_spots(spots_root):
    """检查FISH-QUANT spots文件，确定哪些细胞真正有spots数据"""
    print("=== 检查FISH-QUANT spots文件 ===\n")
    
    spot_files = [f for f in os.listdir(spots_root) if f.endswith('.txt') and '_spots' in f]
    cells_with_spots = {}
    
    for spot_file in spot_files:
        file_path = os.path.join(spots_root, spot_file)
        # 从文件名提取图像名称
        image_name = spot_file.replace('_spots_250622.txt', '').replace('_DIC', '')
        
        # print(f"检查文件: {spot_file}")
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            current_cell = None
            current_cell_has_spots = False
            
            for line in lines:
                line = line.strip()
                
                # 检测新细胞的开始
                if line.startswith('CELL\t'):
                    # 保存之前细胞的状态
                    if current_cell is not None:
                        if image_name not in cells_with_spots:
                            cells_with_spots[image_name] = {}
                        cells_with_spots[image_name][current_cell] = current_cell_has_spots
                    
                    # 重置新细胞
                    current_cell = line.split('\t')[1]
                    current_cell_has_spots = False
                    
                elif line == 'SPOTS':
                    # 找到SPOTS标记，说明当前细胞有spots数据
                    current_cell_has_spots = True
                    
            # 保存最后一个细胞的状态
            if current_cell is not None:
                if image_name not in cells_with_spots:
                    cells_with_spots[image_name] = {}
                cells_with_spots[image_name][current_cell] = current_cell_has_spots
                
        except Exception as e:
            print(f"  错误：无法读取文件 {spot_file}: {e}")
            continue
        
        # 显示这个文件的结果
        if image_name in cells_with_spots:
            cells_info = cells_with_spots[image_name]
            cells_with = [cell for cell, has_spots in cells_info.items() if has_spots]
            cells_without = [cell for cell, has_spots in cells_info.items() if not has_spots]
            
            # print(f"  有spots的细胞: {cells_with}")
            # print(f"  无spots的细胞: {cells_without}")
    
    return cells_with_spots

def filter_valid_cells(df, cells_with_spots_data):
    """过滤掉没有spots数据的细胞"""
    print("=== 过滤无效细胞 ===\n")
    
    initial_count = len(df)
    valid_rows = []
    
    for _, row in df.iterrows():
        image_name = row['image']
        cell_name = row['cell']
        
        # 检查这个细胞是否有spots数据
        if image_name in cells_with_spots_data:
            if cell_name in cells_with_spots_data[image_name]:
                if cells_with_spots_data[image_name][cell_name]:
                    valid_rows.append(row)
                else:
                    print(f"跳过无spots细胞: {image_name} - {cell_name}")
            else:
                print(f"警告：细胞不在spots文件中: {image_name} - {cell_name}")
        else:
            print(f"警告：图像不在spots文件中: {image_name}")
    
    filtered_df = pd.DataFrame(valid_rows)
    final_count = len(filtered_df)
    
    print(f"\n过滤结果:")
    print(f"  原始数据: {initial_count} 行")
    print(f"  过滤后: {final_count} 行")
    print(f"  移除: {initial_count - final_count} 行 ({(initial_count - final_count)/initial_count*100:.1f}%)")
    
    return filtered_df

def analyze_outlier_cells():
    """分析异常outlier细胞的特征和可能原因 - 修正版本，跳过无spots细胞"""
    
    # 读取数据
    base_dir = 'Y333 ATP6 ATP2'
    outlier_file = os.path.join(base_dir, 'interactive_batch_results', 'atp6_outlier_analysis_1.0_fixed.csv')
    summary_file = os.path.join(base_dir, 'interactive_batch_results', 'atp6_all_cells_summary_fixed.csv')
    spots_root = os.path.join(base_dir, 'atp6_spots')
    
    print("=== 异常Outlier细胞分析 (修正版本) ===\n")
    
    # 首先检查哪些细胞有spots数据
    cells_with_spots_data = check_cells_with_spots(spots_root)
    
    # 读取原始数据
    outlier_df = pd.read_csv(outlier_file)
    summary_df = pd.read_csv(summary_file)
    
    print(f"原始数据统计:")
    print(f"  总outlier spots数量: {len(outlier_df)}")
    print(f"  涉及的细胞数量: {len(outlier_df.groupby(['image', 'cell']))}")
    print(f"  总分析细胞数量: {len(summary_df)}")
    
    # 过滤掉没有spots数据的细胞
    outlier_df_filtered = filter_valid_cells(outlier_df, cells_with_spots_data)
    summary_df_filtered = filter_valid_cells(summary_df, cells_with_spots_data)
    
    print(f"\n过滤后数据统计:")
    print(f"  有效outlier spots数量: {len(outlier_df_filtered)}")
    print(f"  有效异常细胞数量: {len(outlier_df_filtered.groupby(['image', 'cell']))}")
    print(f"  有效分析细胞数量: {len(summary_df_filtered)}")
    
    if len(outlier_df_filtered) == 0:
        print("警告：过滤后没有有效的outlier数据!")
        return None, None
    
    # 分析过滤后的outlier距离分布
    print(f"\n有效Outlier距离统计:")
    print(f"最小距离: {outlier_df_filtered['distance_to_skeleton_um'].min():.2f} μm")
    print(f"最大距离: {outlier_df_filtered['distance_to_skeleton_um'].max():.2f} μm")
    print(f"平均距离: {outlier_df_filtered['distance_to_skeleton_um'].mean():.2f} μm")
    print(f"中位数距离: {outlier_df_filtered['distance_to_skeleton_um'].median():.2f} μm")
    
    # 识别真正的极端异常细胞（距离>10μm的）
    extreme_outliers = outlier_df_filtered[outlier_df_filtered['distance_to_skeleton_um'] > 10]
    print(f"\n真正的极端异常spots（>10μm）: {len(extreme_outliers)}")
    print(f"真正的极端异常细胞: {len(extreme_outliers.groupby(['image', 'cell']))}")
    
    if len(extreme_outliers) == 0:
        print("✅ 好消息：过滤后没有极端异常的细胞了!")
    else:
        print("⚠️  仍有极端异常细胞，需要进一步调查")
        
        # 列出仍然异常的细胞
        extreme_cells = extreme_outliers.groupby(['image', 'cell']).agg({
            'distance_to_skeleton_um': ['count', 'mean', 'min', 'max']
        }).round(2)
        extreme_cells.columns = ['outlier_count', 'mean_distance', 'min_distance', 'max_distance']
        
        print("\n仍然异常的细胞:")
        print(extreme_cells.to_string())
    
    # 分析过滤后的Y翻转模式的影响
    print(f"\n=== Y翻转模式分析 (过滤后) ===")
    
    # 正常细胞的Y翻转分布
    normal_cells = summary_df_filtered[summary_df_filtered['mean_distance'] < 1.0]  # 定义正常细胞
    print(f"正常细胞中Y翻转使用率: {normal_cells['y_flip_used'].mean()*100:.1f}%")
    
    # 异常细胞的Y翻转分布
    outlier_cells_filtered = summary_df_filtered[summary_df_filtered['mean_distance'] > 1.0]  # 重新定义异常细胞
    if len(outlier_cells_filtered) > 0:
        print(f"异常细胞中Y翻转使用率: {outlier_cells_filtered['y_flip_used'].mean()*100:.1f}%")
        print(f"\n异常细胞列表:")
        print(outlier_cells_filtered[['mean_distance', 'y_flip_used', 'num_spots']].to_string())
    
    # 分析图像级别的模式
    print(f"\n=== 图像级别分析 (过滤后) ===")
    
    # 按图像统计真实异常细胞数量
    if len(outlier_df_filtered) > 0:
        outlier_by_image = outlier_df_filtered.groupby('image').agg({
            'cell': 'nunique',
            'distance_to_skeleton_um': 'mean'
        }).rename(columns={'cell': 'outlier_cells_count', 'distance_to_skeleton_um': 'avg_outlier_distance'})
    else:
        outlier_by_image = pd.DataFrame()
    
    total_by_image = summary_df_filtered.groupby('image').agg({
        'cell': 'nunique',
        'mean_distance': 'mean'
    }).rename(columns={'cell': 'total_valid_cells', 'mean_distance': 'avg_distance'})
    
    image_analysis = outlier_by_image.join(total_by_image, how='right').fillna(0)
    if len(outlier_by_image) > 0:
        image_analysis['outlier_ratio'] = image_analysis['outlier_cells_count'] / image_analysis['total_valid_cells']
    else:
        image_analysis['outlier_ratio'] = 0
    
    print("每个图像的有效细胞统计:")
    print(image_analysis.round(3).to_string())
    
    # 生成对比可视化
    create_comparison_visualizations(outlier_df, outlier_df_filtered, summary_df, summary_df_filtered, base_dir)
    
    return outlier_df_filtered, summary_df_filtered

def create_comparison_visualizations(outlier_df_orig, outlier_df_filtered, summary_df_orig, summary_df_filtered, base_dir):
    """创建过滤前后的对比可视化"""
    
    output_dir = os.path.join(base_dir, 'interactive_batch_results')
    
    # 1. 距离分布对比 - 过滤前后
    plt.figure(figsize=(15, 10))
    
    # 上排：原始数据
    plt.subplot(2, 3, 1)
    plt.hist(summary_df_orig['mean_distance'], bins=50, alpha=0.7, label='All cells', color='blue')
    plt.hist(summary_df_orig[summary_df_orig['mean_distance'] > 1]['mean_distance'], bins=30, alpha=0.7, label='Outlier cells (>1μm)', color='red')
    plt.xlabel('Mean distance (μm)')
    plt.ylabel('Number of cells')
    plt.title('Distance Distribution (Original)')
    plt.legend()
    plt.yscale('log')
    
    plt.subplot(2, 3, 2)
    extreme_orig = summary_df_orig[summary_df_orig['mean_distance'] > 10]
    if len(extreme_orig) > 0:
        plt.bar(['Normal', 'Extreme'], [len(summary_df_orig) - len(extreme_orig), len(extreme_orig)], color=['blue', 'red'])
        plt.ylabel('Number of cells')
        plt.title(f'Extreme Outliers (Original)\n{len(extreme_orig)} cells')
    else:
        plt.text(0.5, 0.5, 'No extreme outliers', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Extreme Outliers (Original)')
    
    plt.subplot(2, 3, 3)
    if len(outlier_df_orig) > 0:
        plt.hist(outlier_df_orig['distance_to_skeleton_um'], bins=30, alpha=0.7, color='red')
        plt.xlabel('Distance to skeleton (μm)')
        plt.ylabel('Number of spots')
        plt.title('Outlier Spots Distance (Original)')
        plt.yscale('log')
    
    # 下排：过滤后数据
    plt.subplot(2, 3, 4)
    if len(summary_df_filtered) > 0:
        plt.hist(summary_df_filtered['mean_distance'], bins=50, alpha=0.7, label='Valid cells', color='green')
        outlier_filtered = summary_df_filtered[summary_df_filtered['mean_distance'] > 1]
        if len(outlier_filtered) > 0:
            plt.hist(outlier_filtered['mean_distance'], bins=30, alpha=0.7, label='Outlier cells (>1μm)', color='orange')
        plt.xlabel('Mean distance (μm)')
        plt.ylabel('Number of cells')
        plt.title('Distance Distribution (Filtered)')
        plt.legend()
        plt.yscale('log')
    
    plt.subplot(2, 3, 5)
    if len(summary_df_filtered) > 0:
        extreme_filtered = summary_df_filtered[summary_df_filtered['mean_distance'] > 10]
        plt.bar(['Normal', 'Extreme'], [len(summary_df_filtered) - len(extreme_filtered), len(extreme_filtered)], color=['green', 'orange'])
        plt.ylabel('Number of cells')
        plt.title(f'Extreme Outliers (Filtered)\n{len(extreme_filtered)} cells')
    else:
        plt.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Extreme Outliers (Filtered)')
    
    plt.subplot(2, 3, 6)
    if len(outlier_df_filtered) > 0:
        plt.hist(outlier_df_filtered['distance_to_skeleton_um'], bins=30, alpha=0.7, color='orange')
        plt.xlabel('Distance to skeleton (μm)')
        plt.ylabel('Number of spots')
        plt.title('Outlier Spots Distance (Filtered)')
        plt.yscale('log')
    else:
        plt.text(0.5, 0.5, 'No outlier spots', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Outlier Spots Distance (Filtered)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'outlier_analysis_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n可视化图表已保存到: {os.path.join(output_dir, 'outlier_analysis_comparison.png')}")

def check_coordinate_mappings():
    """检查异常细胞的坐标映射文件"""
    
    base_dir = 'Y333 ATP6 ATP2'
    skeleton_root = os.path.join(base_dir, 'extracted_cells')
    
    # 读取outlier数据
    outlier_file = os.path.join(base_dir, 'interactive_batch_results', 'atp6_outlier_analysis_1.0_fixed.csv')
    outlier_df = pd.read_csv(outlier_file)
    
    # 获取极端异常细胞
    extreme_outliers = outlier_df[outlier_df['distance_to_skeleton_um'] > 10]
    extreme_cells = extreme_outliers.groupby(['image', 'cell']).first().index.tolist()
    
    print(f"=== 检查{len(extreme_cells)}个极端异常细胞的坐标映射 ===\n")
    
    for image_name, cell_name in extreme_cells[:5]:  # 检查前5个
        mapping_file = os.path.join(skeleton_root, image_name, 'coordinate_mapping.json')
        
        print(f"{image_name} - {cell_name}:")
        
        if os.path.exists(mapping_file):
            with open(mapping_file, 'r') as f:
                mapping_data = json.load(f)
            
            cell_number = int(cell_name.split("_")[1])
            cell_id_str = f"cell_{cell_number:03d}"
            
            if cell_id_str in mapping_data:
                crop_info = mapping_data[cell_id_str]['crop_region']
                print(f"  Crop region: x_offset={crop_info['x_offset']}, y_offset={crop_info['y_offset']}")
                print(f"               x_start={crop_info['x_start']}, x_end={crop_info['x_end']}")
                print(f"               y_start={crop_info['y_start']}, y_end={crop_info['y_end']}")
                
                # 计算裁剪区域大小
                crop_width = crop_info['x_end'] - crop_info['x_start']
                crop_height = crop_info['y_end'] - crop_info['y_start']
                print(f"  Crop size: {crop_width} x {crop_height} pixels")
                
                # 检查是否有异常的偏移
                if abs(crop_info['x_offset']) > 1000 or abs(crop_info['y_offset']) > 1000:
                    print("  ⚠️  异常大的偏移量!")
                
                if crop_width > 1000 or crop_height > 1000:
                    print("  ⚠️  异常大的裁剪区域!")
                    
            else:
                print(f"  ❌ 未找到{cell_id_str}的映射信息")
        else:
            print("  ❌ 坐标映射文件不存在")
        
        print()

if __name__ == "__main__":
    # 运行修正的分析
    outlier_df_filtered, summary_df_filtered = analyze_outlier_cells()
    
    if outlier_df_filtered is not None:
        print(f"\n=== 总结 ===")
        print("修正分析结果:")
        print("1. ✅ 已跳过所有没有FISH-QUANT spots数据的细胞")
        print("2. ✅ 避免了因为数据错配导致的假性异常细胞")
        print("3. ✅ 现在的分析结果反映了真实的对齐质量")
        
        if len(outlier_df_filtered) == 0:
            print("4. ✅ 太棒了！修正后没有任何outlier细胞")
        else:
            print(f"4. ⚠️  仍有{len(outlier_df_filtered.groupby(['image', 'cell']))}个细胞需要检查")
    else:
        print("⚠️ 无法完成分析，请检查数据文件")
    
    print("\n" + "="*60)
    check_coordinate_mappings()
    
    print(f"\n=== 总结 ===")
    print("异常细胞的可能原因:")
    print("1. Y轴翻转设置错误")
    print("2. 坐标映射文件中的裁剪区域信息错误")
    print("3. 细胞边界检测失败，导致错误的偏移量")
    print("4. 某些图像的预处理步骤出现问题")
    print("5. FISH-QUANT分析结果中的坐标单位或比例错误") 