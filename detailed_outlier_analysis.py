import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats

def analyze_outlier_characteristics():
    """深入分析outlier细胞的特征 - 只使用outlier数据"""
    
    base_dir = 'Y333 ATP6 ATP2'
    outlier_file = os.path.join(base_dir, 'interactive_batch_results', 'atp6_outlier_analysis_1.0_fixed.csv')
    distances_file = os.path.join(base_dir, 'interactive_batch_results', 'atp6_all_cells_distances_fixed.csv')
    
    print("=== 深入分析Outlier细胞特征 (基于距离数据) ===\n")
    
    # 读取数据
    outlier_df = pd.read_csv(outlier_file)
    
    # 尝试读取距离数据作为总体统计的参考
    distances_df = None
    if os.path.exists(distances_file):
        try:
            distances_df = pd.read_csv(distances_file)
            # 检查是否有必要的列用于匹配
            if 'image' in distances_df.columns and 'cell' in distances_df.columns:
                print("已找到完整的距离数据文件，将用于总体统计")
            else:
                print("距离数据文件格式不完整，只有距离值，将用于基本统计对比")
        except Exception as e:
            print(f"读取距离数据文件失败: {e}")
            distances_df = None
    else:
        print("未找到距离数据文件，将只分析outlier数据")
    
    # 移除问题图像数据
    problem_image = 'yWL333_cy3_ATP2_cy5_ATP6MS2_2_s1'
    outlier_df = outlier_df[outlier_df['image'] != problem_image]
    # 只有当distances_df有image列时才过滤
    if distances_df is not None and 'image' in distances_df.columns:
        distances_df = distances_df[distances_df['image'] != problem_image]
    
    print(f"分析数据统计:")
    print(f"  有效outlier spots: {len(outlier_df)}")
    print(f"  有效异常细胞: {len(outlier_df.groupby(['image', 'cell']))}")
    if distances_df is not None:
        print(f"  总分析spots: {len(distances_df)}")
        if 'image' in distances_df.columns and 'cell' in distances_df.columns:
            print(f"  总分析细胞: {len(distances_df.groupby(['image', 'cell']))}")
        else:
            print(f"  (distances数据无细胞分组信息)")
    
    # 1. 分析outlier距离分布
    analyze_outlier_distance_distribution(outlier_df, distances_df)
    
    # 2. 分析outlier细胞的特征
    analyze_outlier_cell_characteristics(outlier_df)
    
    # 3. 分析图像级别的outlier模式
    analyze_outlier_image_patterns(outlier_df)
    
    # 4. 识别需要检查的具体细胞
    identify_outlier_cells_for_inspection(outlier_df)
    
    # 5. 生成可视化
    create_outlier_visualizations(outlier_df, distances_df, base_dir)
    
    return outlier_df, distances_df

def analyze_outlier_distance_distribution(outlier_df, distances_df):
    """分析outlier距离分布"""
    
    print(f"\n=== 1. Outlier距离分布分析 ===")
    
    # 分析outlier距离
    outlier_distances = outlier_df['distance_to_skeleton_um']
    
    print(f"Outlier spots距离统计:")
    print(f"  总outlier spots: {len(outlier_df)}")
    print(f"  平均距离: {outlier_distances.mean():.3f} μm")
    print(f"  中位数: {outlier_distances.median():.3f} μm")
    print(f"  标准差: {outlier_distances.std():.3f} μm")
    print(f"  最小距离: {outlier_distances.min():.3f} μm")
    print(f"  最大距离: {outlier_distances.max():.3f} μm")
    print(f"  95%分位数: {outlier_distances.quantile(0.95):.3f} μm")
    print(f"  99%分位数: {outlier_distances.quantile(0.99):.3f} μm")
    
    # 距离分布
    distance_ranges = [
        (1.0, 1.2),
        (1.2, 1.5),
        (1.5, 2.0),
        (2.0, 3.0),
        (3.0, float('inf'))
    ]
    
    print(f"\nOutlier距离分布:")
    for min_dist, max_dist in distance_ranges:
        if max_dist == float('inf'):
            count = len(outlier_df[outlier_df['distance_to_skeleton_um'] >= min_dist])
            label = f"≥{min_dist:.1f}μm"
        else:
            count = len(outlier_df[(outlier_df['distance_to_skeleton_um'] >= min_dist) & 
                                  (outlier_df['distance_to_skeleton_um'] < max_dist)])
            label = f"{min_dist:.1f}-{max_dist:.1f}μm"
        
        percentage = count / len(outlier_df) * 100
        print(f"  {label}: {count} spots ({percentage:.1f}%)")
    
    # 如果有总体距离数据，进行对比
    if distances_df is not None:
        print(f"\n总体数据对比:")
        # 检查距离列名
        if 'distance_to_skeleton_um' in distances_df.columns:
            distance_col = 'distance_to_skeleton_um'
        elif 'distance' in distances_df.columns:
            distance_col = 'distance'
        else:
            distance_col = distances_df.columns[0]  # 使用第一列
            
        all_distances = distances_df[distance_col]
        print(f"  总spots数: {len(distances_df)}")
        print(f"  总体平均距离: {all_distances.mean():.3f} μm")
        print(f"  总体中位数: {all_distances.median():.3f} μm")
        print(f"  outlier比例: {len(outlier_df)/len(distances_df)*100:.2f}%")

def analyze_outlier_cell_characteristics(outlier_df):
    """分析outlier细胞特征"""
    
    print(f"\n=== 2. Outlier细胞特征分析 ===")
    
    # 按细胞分组分析
    cell_stats = outlier_df.groupby(['image', 'cell']).agg({
        'distance_to_skeleton_um': ['count', 'mean', 'min', 'max', 'std']
    }).round(3)
    cell_stats.columns = ['outlier_spots_count', 'mean_outlier_distance', 'min_distance', 'max_distance', 'distance_std']
    cell_stats = cell_stats.fillna(0)
    
    print(f"异常细胞总数: {len(cell_stats)} 个")
    print(f"总outlier spots: {len(outlier_df)} 个")
    print(f"平均每个异常细胞的outlier spots数: {len(outlier_df)/len(cell_stats):.1f}")
    
    # 按每个细胞的outlier spots数量分析
    spots_counts = cell_stats['outlier_spots_count']
    print(f"\n每个异常细胞的outlier spots数量分布:")
    print(f"  平均: {spots_counts.mean():.1f}")
    print(f"  中位数: {spots_counts.median():.1f}")
    print(f"  标准差: {spots_counts.std():.1f}")
    print(f"  最小: {spots_counts.min()}")
    print(f"  最大: {spots_counts.max()}")
    
    # 按outlier spots数量分类
    print(f"\n按outlier spots数量分类:")
    for i in range(1, spots_counts.max() + 1):
        count = len(cell_stats[cell_stats['outlier_spots_count'] == i])
        percentage = count / len(cell_stats) * 100
        print(f"  {i}个outlier spots: {count} 细胞 ({percentage:.1f}%)")
    
    # 按平均距离分析
    mean_distances = cell_stats['mean_outlier_distance']
    print(f"\n各异常细胞的平均outlier距离:")
    print(f"  平均: {mean_distances.mean():.3f} μm")
    print(f"  中位数: {mean_distances.median():.3f} μm")
    print(f"  标准差: {mean_distances.std():.3f} μm")
    print(f"  最小: {mean_distances.min():.3f} μm")
    print(f"  最大: {mean_distances.max():.3f} μm")
    
    # 识别最严重的异常细胞
    severe_cells = cell_stats[cell_stats['mean_outlier_distance'] > 1.5]
    print(f"\n严重异常细胞 (平均距离>1.5μm): {len(severe_cells)} 个")
    if len(severe_cells) > 0:
        print(f"  占比: {len(severe_cells)/len(cell_stats)*100:.1f}%")
        print("  具体细胞:")
        for idx, (image, cell) in enumerate(severe_cells.index):
            stats_row = severe_cells.loc[(image, cell)]
            print(f"    {image} - {cell}: 平均{stats_row['mean_outlier_distance']:.3f}μm, {stats_row['outlier_spots_count']:.0f}个outlier spots")
    
    return cell_stats

def analyze_outlier_image_patterns(outlier_df):
    """分析图像级别的outlier模式"""
    
    print(f"\n=== 3. 图像级别Outlier模式分析 ===")
    
    # 按图像统计outlier情况
    image_outlier_stats = outlier_df.groupby('image').agg({
        'distance_to_skeleton_um': ['count', 'mean', 'min', 'max', 'std'],
        'cell': 'nunique'
    }).round(3)
    image_outlier_stats.columns = ['total_outlier_spots', 'avg_outlier_distance', 'min_distance', 'max_distance', 'distance_std', 'outlier_cells_count']
    
    # 按outlier spots数量排序
    image_outlier_stats = image_outlier_stats.sort_values('total_outlier_spots', ascending=False)
    
    print(f"按outlier spots数量排序的图像:")
    print(f"图像数量: {len(image_outlier_stats)}")
    print(f"总outlier spots: {image_outlier_stats['total_outlier_spots'].sum()}")
    print(f"总异常细胞: {image_outlier_stats['outlier_cells_count'].sum()}")
    print()
    
    # 显示前10个最严重的图像
    print("前10个outlier最多的图像:")
    top_images = image_outlier_stats.head(10)
    for image in top_images.index:
        stats = top_images.loc[image]
        print(f"  {image}:")
        print(f"    Outlier spots: {stats['total_outlier_spots']:.0f}")
        print(f"    异常细胞: {stats['outlier_cells_count']:.0f}")
        print(f"    平均距离: {stats['avg_outlier_distance']:.3f}μm")
        print(f"    距离范围: {stats['min_distance']:.3f}-{stats['max_distance']:.3f}μm")
        print()
    
    # 分析图像系列模式
    print("按图像系列分析:")
    series_patterns = {}
    for image in image_outlier_stats.index:
        # 提取系列信息 (例如: yWL333_cy3_ATP2_cy5_ATP6MS2_1_s1 -> 1)
        try:
            series = image.split('_')[-2]  # 获取倒数第二个部分，应该是系列号
            if series not in series_patterns:
                series_patterns[series] = {
                    'images': [],
                    'total_outlier_spots': 0,
                    'total_outlier_cells': 0
                }
            series_patterns[series]['images'].append(image)
            series_patterns[series]['total_outlier_spots'] += image_outlier_stats.loc[image, 'total_outlier_spots']
            series_patterns[series]['total_outlier_cells'] += image_outlier_stats.loc[image, 'outlier_cells_count']
        except:
            continue
    
    for series, data in sorted(series_patterns.items()):
        avg_spots_per_image = data['total_outlier_spots'] / len(data['images'])
        avg_cells_per_image = data['total_outlier_cells'] / len(data['images'])
        print(f"  系列 {series}: {len(data['images'])} 图像")
        print(f"    总outlier spots: {data['total_outlier_spots']}")
        print(f"    总异常细胞: {data['total_outlier_cells']}")
        print(f"    平均每图像: {avg_spots_per_image:.1f} spots, {avg_cells_per_image:.1f} 细胞")
        print()
    
    return image_outlier_stats

def identify_outlier_cells_for_inspection(outlier_df):
    """识别需要人工检查的具体outlier细胞"""
    
    print(f"\n=== 4. 需要检查的具体Outlier细胞 ===")
    
    # 先按细胞分组获取统计信息
    cell_stats = outlier_df.groupby(['image', 'cell']).agg({
        'distance_to_skeleton_um': ['count', 'mean', 'min', 'max']
    }).round(3)
    cell_stats.columns = ['outlier_spots_count', 'mean_outlier_distance', 'min_distance', 'max_distance']
    
    print(f"异常细胞总数: {len(cell_stats)}")
    
    # 按优先级分类
    print(f"\n按优先级分类:")
    
    # 高优先级：平均距离>2.0μm或最大距离>3.0μm
    high_priority = cell_stats[(cell_stats['mean_outlier_distance'] > 2.0) | 
                              (cell_stats['max_distance'] > 3.0)]
    print(f"1. 高优先级: {len(high_priority)} 个细胞")
    print("   (平均距离>2.0μm 或 最大距离>3.0μm)")
    if len(high_priority) > 0:
        print("   具体细胞:")
        for (image, cell), stats in high_priority.iterrows():
            print(f"     {image} - {cell}:")
            print(f"       平均距离: {stats['mean_outlier_distance']:.3f}μm")
            print(f"       距离范围: {stats['min_distance']:.3f}-{stats['max_distance']:.3f}μm")
            print(f"       Outlier spots: {stats['outlier_spots_count']:.0f}个")
            print()
    
    # 中优先级：平均距离1.5-2.0μm
    med_priority = cell_stats[(cell_stats['mean_outlier_distance'] > 1.5) & 
                             (cell_stats['mean_outlier_distance'] <= 2.0) &
                             (cell_stats['max_distance'] <= 3.0)]
    print(f"2. 中优先级: {len(med_priority)} 个细胞")
    print("   (平均距离1.5-2.0μm)")
    if len(med_priority) > 0:
        print("   具体细胞:")
        for (image, cell), stats in med_priority.iterrows():
            print(f"     {image} - {cell}: 平均{stats['mean_outlier_distance']:.3f}μm, {stats['outlier_spots_count']:.0f}个outlier spots")
    
    # 低优先级：平均距离1.0-1.5μm
    low_priority = cell_stats[(cell_stats['mean_outlier_distance'] > 1.0) & 
                             (cell_stats['mean_outlier_distance'] <= 1.5)]
    print(f"\n3. 低优先级: {len(low_priority)} 个细胞")
    print("   (平均距离1.0-1.5μm)")
    
    # 特殊关注：只有1个outlier spot的细胞
    single_spot_outliers = cell_stats[cell_stats['outlier_spots_count'] == 1]
    print(f"\n4. 特殊关注 (只有1个outlier spot): {len(single_spot_outliers)} 个细胞")
    print("   这些细胞可能是假阳性，建议优先检查:")
    if len(single_spot_outliers) > 0:
        for (image, cell), stats in single_spot_outliers.iterrows():
            print(f"     {image} - {cell}: 距离{stats['mean_outlier_distance']:.3f}μm")
    
    # 多outlier spots的严重异常细胞
    multiple_spots_outliers = cell_stats[cell_stats['outlier_spots_count'] >= 3]
    print(f"\n5. 多outlier spots细胞 (≥3个): {len(multiple_spots_outliers)} 个细胞")
    print("   这些细胞的异常模式可能更可靠:")
    if len(multiple_spots_outliers) > 0:
        for (image, cell), stats in multiple_spots_outliers.iterrows():
            print(f"     {image} - {cell}: 平均{stats['mean_outlier_distance']:.3f}μm, {stats['outlier_spots_count']:.0f}个outlier spots")
    
    return cell_stats

def create_outlier_visualizations(outlier_df, distances_df, base_dir):
    """创建基于outlier数据的可视化图表"""
    
    output_dir = os.path.join(base_dir, 'interactive_batch_results')
    
    plt.figure(figsize=(16, 12))
    
    # 1. Outlier距离分布直方图
    plt.subplot(3, 4, 1)
    plt.hist(outlier_df['distance_to_skeleton_um'], bins=30, alpha=0.7, color='red', label='Outlier spots')
    plt.axvline(outlier_df['distance_to_skeleton_um'].mean(), color='darkred', linestyle='--', label='Mean')
    plt.axvline(outlier_df['distance_to_skeleton_um'].median(), color='orange', linestyle='--', label='Median')
    plt.xlabel('Distance (μm)')
    plt.ylabel('Number of outlier spots')
    plt.title('Outlier Distance Distribution')
    plt.legend()
    
    # 2. 总体距离分布对比（如果有数据）
    plt.subplot(3, 4, 2)
    if distances_df is not None:
        # 检查distances_df的列结构
        if 'distance_to_skeleton_um' in distances_df.columns:
            distance_col = 'distance_to_skeleton_um'
        elif 'distance' in distances_df.columns:
            distance_col = 'distance'
        else:
            distance_col = distances_df.columns[0]  # 使用第一列
        
        plt.hist(distances_df[distance_col], bins=100, alpha=0.5, color='blue', label='All spots', density=True)
        plt.hist(outlier_df['distance_to_skeleton_um'], bins=30, alpha=0.7, color='red', label='Outlier spots', density=True)
        plt.axvline(1.0, color='black', linestyle='--', label='Threshold')
        plt.xlabel('Distance (μm)')
        plt.ylabel('Density')
        plt.title('Distance Distribution Comparison')
        plt.legend()
        plt.yscale('log')
    else:
        plt.text(0.5, 0.5, 'No overall\ndistance data', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Overall Distance Distribution')
    
    # 3. 按细胞的outlier spots数量
    plt.subplot(3, 4, 3)
    cell_spot_counts = outlier_df.groupby(['image', 'cell']).size()
    plt.hist(cell_spot_counts, bins=range(1, cell_spot_counts.max()+2), alpha=0.7, color='orange')
    plt.xlabel('Outlier spots per cell')
    plt.ylabel('Number of cells')
    plt.title('Outlier Spots per Cell')
    
    # 4. 按图像的outlier分布
    plt.subplot(3, 4, 4)
    image_spot_counts = outlier_df.groupby('image').size()
    plt.bar(range(len(image_spot_counts)), sorted(image_spot_counts.values, reverse=True), color='coral')
    plt.xlabel('Image (ranked by outlier count)')
    plt.ylabel('Number of outlier spots')
    plt.title('Outlier Spots by Image')
    
    # 5. 距离vs细胞散点图
    plt.subplot(3, 4, 5)
    cell_stats = outlier_df.groupby(['image', 'cell']).agg({
        'distance_to_skeleton_um': ['count', 'mean']
    })
    cell_stats.columns = ['spot_count', 'mean_distance']
    
    plt.scatter(cell_stats['spot_count'], cell_stats['mean_distance'], alpha=0.7, color='red')
    plt.xlabel('Number of outlier spots')
    plt.ylabel('Mean outlier distance (μm)')
    plt.title('Cell Outlier Characteristics')
    
    # 6. 累积分布图
    plt.subplot(3, 4, 6)
    sorted_distances = np.sort(outlier_df['distance_to_skeleton_um'])
    cumulative = np.arange(1, len(sorted_distances) + 1) / len(sorted_distances)
    plt.plot(sorted_distances, cumulative, 'r-', linewidth=2, label='Outlier spots')
    plt.axvline(sorted_distances[int(0.5 * len(sorted_distances))], color='orange', linestyle='--', label='Median')
    plt.axvline(sorted_distances[int(0.95 * len(sorted_distances))], color='red', linestyle='--', label='95th percentile')
    plt.xlabel('Distance (μm)')
    plt.ylabel('Cumulative probability')
    plt.title('Outlier Distance CDF')
    plt.legend()
    
    # 7. 按图像系列的outlier分布
    plt.subplot(3, 4, 7)
    series_data = {}
    for image in outlier_df['image'].unique():
        try:
            series = image.split('_')[-2]  # 获取系列号
            if series not in series_data:
                series_data[series] = 0
            series_data[series] += len(outlier_df[outlier_df['image'] == image])
        except:
            continue
    
    if series_data:
        plt.bar(series_data.keys(), series_data.values(), color='lightcoral')
        plt.xlabel('Series')
        plt.ylabel('Total outlier spots')
        plt.title('Outlier Spots by Series')
    
    # 8. 距离范围分布
    plt.subplot(3, 4, 8)
    distance_ranges = [(1.0, 1.2), (1.2, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, float('inf'))]
    range_counts = []
    range_labels = []
    
    for min_dist, max_dist in distance_ranges:
        if max_dist == float('inf'):
            count = len(outlier_df[outlier_df['distance_to_skeleton_um'] >= min_dist])
            label = f'≥{min_dist:.1f}μm'
        else:
            count = len(outlier_df[(outlier_df['distance_to_skeleton_um'] >= min_dist) & 
                                  (outlier_df['distance_to_skeleton_um'] < max_dist)])
            label = f'{min_dist:.1f}-{max_dist:.1f}μm'
        range_counts.append(count)
        range_labels.append(label)
    
    plt.bar(range_labels, range_counts, color='salmon')
    plt.xlabel('Distance range')
    plt.ylabel('Number of spots')
    plt.title('Outlier Distance Ranges')
    plt.xticks(rotation=45)
    
    # 9. 箱线图：每个细胞的outlier距离分布
    plt.subplot(3, 4, 9)
    cell_distances = []
    cell_labels = []
    
    # 只显示前10个有最多outlier spots的细胞
    top_cells = outlier_df.groupby(['image', 'cell']).size().nlargest(10)
    for (image, cell), count in top_cells.items():
        cell_data = outlier_df[(outlier_df['image'] == image) & (outlier_df['cell'] == cell)]
        cell_distances.append(cell_data['distance_to_skeleton_um'].values)
        cell_labels.append(f'{cell}({count})')
    
    if cell_distances:
        plt.boxplot(cell_distances, labels=cell_labels)
        plt.xlabel('Cell (outlier count)')
        plt.ylabel('Distance (μm)')
        plt.title('Top Outlier Cells Distribution')
        plt.xticks(rotation=45)
    
    # 10. 优先级分布饼图
    plt.subplot(3, 4, 10)
    cell_mean_distances = outlier_df.groupby(['image', 'cell'])['distance_to_skeleton_um'].mean()
    
    high_priority = len(cell_mean_distances[cell_mean_distances > 2.0])
    med_priority = len(cell_mean_distances[(cell_mean_distances > 1.5) & (cell_mean_distances <= 2.0)])
    low_priority = len(cell_mean_distances[(cell_mean_distances > 1.0) & (cell_mean_distances <= 1.5)])
    
    if high_priority + med_priority + low_priority > 0:
        plt.pie([high_priority, med_priority, low_priority], 
                labels=[f'High (>2.0μm)\n{high_priority}', f'Med (1.5-2.0μm)\n{med_priority}', f'Low (1.0-1.5μm)\n{low_priority}'],
                autopct='%1.1f%%', colors=['red', 'orange', 'yellow'])
        plt.title('Outlier Cell Priority')
    
    # 11-12. 预留空间
    for i in range(11, 13):
        plt.subplot(3, 4, i)
        plt.text(0.5, 0.5, f'Reserved\nfor future\nanalysis', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title(f'Analysis {i}')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'outlier_focused_analysis.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nOutlier专项可视化图表已保存到: {os.path.join(output_dir, 'outlier_focused_analysis.png')}")

if __name__ == "__main__":
    outlier_df, distances_df = analyze_outlier_characteristics()
    
    print(f"\n=== 最终建议 ===")
    print("1. 🎯 基于outlier数据的分析结果:")
    print("   - 重点关注平均距离>2.0μm的细胞")
    print("   - 特别检查只有1个outlier spot的细胞（可能假阳性）")
    print("   - 多outlier spots的细胞可能需要进一步调查")
    
    print("2. 🔍 优先检查顺序:")
    print("   - 高优先级：平均距离>2.0μm或最大距离>3.0μm")
    print("   - 中优先级：平均距离1.5-2.0μm")
    print("   - 特殊关注：只有1个outlier spot的细胞")
    
    print("3. 📊 需要查看原始图像的情况:")
    print("   - 验证高优先级细胞的FISH信号质量")
    print("   - 检查单spot异常细胞是否为假阳性")
    print("   - 确认多outlier spots细胞的对齐问题")
    print("   - 查看outlier最多的图像是否有系统性问题") 