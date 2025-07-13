import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

def load_outlier_data(file_path):
    """Load and validate the outlier CSV file"""
    print("=== Loading Outlier Data ===\n")
    
    if not os.path.exists(file_path):
        print(f"❌ Error: File not found at {file_path}")
        return None
    
    try:
        df = pd.read_csv(file_path)
        print(f"✅ Successfully loaded data from: {os.path.basename(file_path)}")
        print(f"   Total outlier spots: {len(df)}")
        print(f"   Columns: {list(df.columns)}")
        
        # Basic validation
        required_cols = ['image', 'cell', 'distance_to_skeleton_um', 'threshold_um']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"⚠️  Warning: Missing required columns: {missing_cols}")
        
        return df
    
    except Exception as e:
        print(f"❌ Error loading file: {e}")
        return None

def analyze_basic_statistics(df):
    """Analyze basic statistics of the outlier data"""
    print("\n=== Basic Statistics ===\n")
    
    # Overall statistics
    print("📊 Overall Outlier Statistics:")
    print(f"   Total outlier spots: {len(df)}")
    print(f"   Unique cells with outliers: {df.groupby(['image', 'cell']).ngroups}")
    print(f"   Unique images: {df['image'].nunique()}")
    print(f"   Threshold used: {df['threshold_um'].iloc[0]} μm")
    
    # Distance statistics
    distances = df['distance_to_skeleton_um']
    print(f"\n📏 Distance Statistics:")
    print(f"   Min distance: {distances.min():.3f} μm")
    print(f"   Max distance: {distances.max():.3f} μm")
    print(f"   Mean distance: {distances.mean():.3f} μm")
    print(f"   Median distance: {distances.median():.3f} μm")
    print(f"   Std deviation: {distances.std():.3f} μm")
    
    # Percentiles
    print(f"\n📈 Distance Percentiles:")
    for p in [75, 90, 95, 99]:
        print(f"   {p}th percentile: {np.percentile(distances, p):.3f} μm")
    
    # Extreme outliers
    extreme_spots = df[df['distance_to_skeleton_um'] > 2.0]
    if len(extreme_spots) > 0:
        print(f"\n🚨 Extreme Outliers (>2.0 μm): {len(extreme_spots)} spots")
        extreme_cells = extreme_spots.groupby(['image', 'cell']).ngroups
        print(f"   Affecting {extreme_cells} cells")
    else:
        print(f"\n✅ No extreme outliers (>2.0 μm) found")

def analyze_by_cell(df):
    """Analyze outliers grouped by cell"""
    print("\n=== Cell-Level Analysis ===\n")
    
    # Group by cell
    cell_stats = df.groupby(['image', 'cell']).agg({
        'distance_to_skeleton_um': ['count', 'mean', 'max', 'std'],
        'y_flip_used': 'first',
        'z_flip_used': 'first'
    }).round(3)
    
    # Flatten column names
    cell_stats.columns = ['spot_count', 'mean_distance', 'max_distance', 'std_distance', 'y_flip', 'z_flip']
    cell_stats = cell_stats.sort_values('mean_distance', ascending=False)
    
    print("🔍 Top 10 Most Problematic Cells:")
    print(cell_stats.head(10).to_string())
    
    # Statistics about cells
    print(f"\n📊 Cell Distribution:")
    print(f"   Cells with 1 outlier spot: {sum(cell_stats['spot_count'] == 1)}")
    print(f"   Cells with 2-5 outlier spots: {sum((cell_stats['spot_count'] >= 2) & (cell_stats['spot_count'] <= 5))}")
    print(f"   Cells with >5 outlier spots: {sum(cell_stats['spot_count'] > 5)}")
    
    # Y-flip analysis
    if 'y_flip' in cell_stats.columns:
        y_flip_usage = cell_stats['y_flip'].value_counts()
        print(f"\n🔄 Y-flip Usage in Outlier Cells:")
        for flip_state, count in y_flip_usage.items():
            print(f"   Y-flip {flip_state}: {count} cells ({count/len(cell_stats)*100:.1f}%)")
    
    return cell_stats

def analyze_by_image(df):
    """Analyze outliers grouped by image"""
    print("\n=== Image-Level Analysis ===\n")
    
    # Group by image
    image_stats = df.groupby('image').agg({
        'cell': 'nunique',
        'distance_to_skeleton_um': ['count', 'mean', 'max'],
        'y_flip_used': lambda x: x.sum() / len(x) * 100  # percentage of spots with y_flip
    }).round(3)
    
    # Flatten column names
    image_stats.columns = ['unique_cells', 'total_spots', 'mean_distance', 'max_distance', 'y_flip_percent']
    image_stats = image_stats.sort_values('mean_distance', ascending=False)
    
    print("🖼️  Images Ranked by Average Outlier Distance:")
    print(image_stats.to_string())
    
    # Find most problematic images
    worst_images = image_stats.head(3)
    print(f"\n🚨 Top 3 Most Problematic Images:")
    for idx, (image_name, stats) in enumerate(worst_images.iterrows(), 1):
        print(f"   {idx}. {image_name}")
        print(f"      - {stats['unique_cells']} cells with outliers")
        print(f"      - {stats['total_spots']} total outlier spots")
        print(f"      - Average distance: {stats['mean_distance']:.3f} μm")
    
    return image_stats

def create_visualizations(df, cell_stats, image_stats, output_dir):
    """Create comprehensive visualizations"""
    print(f"\n=== Creating Visualizations ===\n")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create a comprehensive figure
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Distance distribution histogram
    plt.subplot(3, 4, 1)
    plt.hist(df['distance_to_skeleton_um'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(df['distance_to_skeleton_um'].mean(), color='red', linestyle='--', label=f'Mean: {df["distance_to_skeleton_um"].mean():.3f}μm')
    plt.xlabel('Distance to Skeleton (μm)')
    plt.ylabel('Number of Spots')
    plt.title('Distribution of Outlier Distances')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Box plot by image (top 10)
    plt.subplot(3, 4, 2)
    top_images = image_stats.head(10).index
    plot_data = df[df['image'].isin(top_images)]
    sns.boxplot(data=plot_data, x='distance_to_skeleton_um', y='image', orient='h')
    plt.xlabel('Distance to Skeleton (μm)')
    plt.title('Distance Distribution by Image (Top 10)')
    plt.tight_layout()
    
    # 3. Spots per cell histogram
    plt.subplot(3, 4, 3)
    plt.hist(cell_stats['spot_count'], bins=range(1, cell_stats['spot_count'].max()+2), 
             alpha=0.7, color='lightgreen', edgecolor='black')
    plt.xlabel('Number of Outlier Spots per Cell')
    plt.ylabel('Number of Cells')
    plt.title('Distribution of Outlier Spots per Cell')
    plt.grid(True, alpha=0.3)
    
    # 4. Y-flip usage
    plt.subplot(3, 4, 4)
    if 'y_flip' in cell_stats.columns:
        y_flip_counts = cell_stats['y_flip'].value_counts()
        plt.pie(y_flip_counts.values, labels=[f'Y-flip {bool(x)}' for x in y_flip_counts.index], 
                autopct='%1.1f%%', startangle=90)
        plt.title('Y-flip Usage in Outlier Cells')
    
    # 5. Scatter plot: Mean distance vs spot count per cell
    plt.subplot(3, 4, 5)
    scatter = plt.scatter(cell_stats['spot_count'], cell_stats['mean_distance'], 
                         alpha=0.6, c=cell_stats['max_distance'], cmap='viridis', s=50)
    plt.colorbar(scatter, label='Max Distance (μm)')
    plt.xlabel('Number of Outlier Spots')
    plt.ylabel('Mean Distance (μm)')
    plt.title('Cell Outlier Profile')
    plt.grid(True, alpha=0.3)
    
    # 6. Distance vs Z coordinate
    plt.subplot(3, 4, 6)
    if 'spot_z_um' in df.columns:
        plt.scatter(df['spot_z_um'], df['distance_to_skeleton_um'], alpha=0.6, s=20)
        plt.xlabel('Z Coordinate (μm)')
        plt.ylabel('Distance to Skeleton (μm)')
        plt.title('Distance vs Z Position')
        plt.grid(True, alpha=0.3)
    
    # 7. Cumulative distribution
    plt.subplot(3, 4, 7)
    sorted_distances = np.sort(df['distance_to_skeleton_um'])
    cumulative_prob = np.arange(1, len(sorted_distances) + 1) / len(sorted_distances)
    plt.plot(sorted_distances, cumulative_prob, linewidth=2)
    plt.axvline(df['threshold_um'].iloc[0], color='red', linestyle='--', label=f'Threshold: {df["threshold_um"].iloc[0]}μm')
    plt.xlabel('Distance to Skeleton (μm)')
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Distribution of Distances')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 8. Image outlier summary
    plt.subplot(3, 4, 8)
    top_5_images = image_stats.head(5)
    plt.barh(range(len(top_5_images)), top_5_images['total_spots'], color='coral')
    plt.yticks(range(len(top_5_images)), [name.split('_')[-1] for name in top_5_images.index])
    plt.xlabel('Total Outlier Spots')
    plt.title('Top 5 Images by Outlier Count')
    plt.grid(True, alpha=0.3)
    
    # 9. Distance vs X,Y coordinates
    plt.subplot(3, 4, 9)
    if 'spot_x_um' in df.columns and 'spot_y_um' in df.columns:
        scatter = plt.scatter(df['spot_x_um'], df['spot_y_um'], 
                             c=df['distance_to_skeleton_um'], cmap='plasma', alpha=0.6, s=20)
        plt.colorbar(scatter, label='Distance (μm)')
        plt.xlabel('X Coordinate (μm)')
        plt.ylabel('Y Coordinate (μm)')
        plt.title('Spatial Distribution of Outliers')
    
    # 10. Distance range per cell
    plt.subplot(3, 4, 10)
    cell_ranges = cell_stats['max_distance'] - cell_stats['mean_distance']
    plt.hist(cell_ranges, bins=20, alpha=0.7, color='orange', edgecolor='black')
    plt.xlabel('Distance Range per Cell (μm)')
    plt.ylabel('Number of Cells')
    plt.title('Variability of Distances within Cells')
    plt.grid(True, alpha=0.3)
    
    # 11. Timeline/sequence analysis if spot_index exists
    plt.subplot(3, 4, 11)
    if 'spot_index' in df.columns:
        plt.scatter(df['spot_index'], df['distance_to_skeleton_um'], alpha=0.5, s=20)
        plt.xlabel('Spot Index')
        plt.ylabel('Distance to Skeleton (μm)')
        plt.title('Distance vs Spot Index')
        plt.grid(True, alpha=0.3)
    
    # 12. Summary statistics text
    plt.subplot(3, 4, 12)
    plt.axis('off')
    summary_text = f"""
SUMMARY STATISTICS

Total Outlier Spots: {len(df)}
Affected Cells: {df.groupby(['image', 'cell']).ngroups}
Affected Images: {df['image'].nunique()}

Distance Stats:
  Mean: {df['distance_to_skeleton_um'].mean():.3f} μm
  Median: {df['distance_to_skeleton_um'].median():.3f} μm
  Max: {df['distance_to_skeleton_um'].max():.3f} μm

Extreme Outliers (>2μm):
  {len(df[df['distance_to_skeleton_um'] > 2.0])} spots
  {df[df['distance_to_skeleton_um'] > 2.0].groupby(['image', 'cell']).ngroups} cells

Most Problematic Cell:
  {cell_stats.index[0][1]} in {cell_stats.index[0][0].split('_')[-1]}
  {cell_stats.iloc[0]['spot_count']} spots, {cell_stats.iloc[0]['mean_distance']:.3f}μm avg
"""
    plt.text(0.1, 0.9, summary_text, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    # Save the comprehensive plot
    output_file = os.path.join(output_dir, 'comprehensive_outlier_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Comprehensive visualization saved: {output_file}")
    
    # Create a focused summary plot
    create_summary_plot(df, cell_stats, output_dir)

def create_summary_plot(df, cell_stats, output_dir):
    """Create a focused summary plot for quick overview"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Distance distribution
    ax1.hist(df['distance_to_skeleton_um'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.axvline(df['distance_to_skeleton_um'].mean(), color='red', linestyle='--', 
                label=f'Mean: {df["distance_to_skeleton_um"].mean():.3f}μm')
    ax1.set_xlabel('Distance to Skeleton (μm)')
    ax1.set_ylabel('Number of Spots')
    ax1.set_title('Distribution of Outlier Distances')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Spots per cell
    ax2.hist(cell_stats['spot_count'], bins=range(1, cell_stats['spot_count'].max()+2), 
             alpha=0.7, color='lightgreen', edgecolor='black')
    ax2.set_xlabel('Outlier Spots per Cell')
    ax2.set_ylabel('Number of Cells')
    ax2.set_title('Outlier Spots Distribution per Cell')
    ax2.grid(True, alpha=0.3)
    
    # 3. Top problematic cells
    top_cells = cell_stats.head(10)
    cell_labels = [f"{idx[1]}\n({idx[0].split('_')[-1]})" for idx in top_cells.index]
    bars = ax3.barh(range(len(top_cells)), top_cells['mean_distance'], color='coral')
    ax3.set_yticks(range(len(top_cells)))
    ax3.set_yticklabels(cell_labels, fontsize=8)
    ax3.set_xlabel('Mean Distance (μm)')
    ax3.set_title('Top 10 Most Problematic Cells')
    ax3.grid(True, alpha=0.3)
    
    # Add values on bars
    for i, bar in enumerate(bars):
        width = bar.get_width()
        ax3.text(width + 0.01, bar.get_y() + bar.get_height()/2, 
                f'{width:.3f}', ha='left', va='center', fontsize=8)
    
    # 4. Summary statistics
    ax4.axis('off')
    stats_text = f"""
OUTLIER ANALYSIS SUMMARY

📊 Data Overview:
   • Total outlier spots: {len(df):,}
   • Affected cells: {df.groupby(['image', 'cell']).ngroups:,}
   • Images with outliers: {df['image'].nunique():,}
   • Threshold used: {df['threshold_um'].iloc[0]} μm

📏 Distance Statistics:
   • Mean: {df['distance_to_skeleton_um'].mean():.3f} μm
   • Median: {df['distance_to_skeleton_um'].median():.3f} μm
   • Max: {df['distance_to_skeleton_um'].max():.3f} μm
   • 95th percentile: {np.percentile(df['distance_to_skeleton_um'], 95):.3f} μm

🚨 Severity Assessment:
   • Extreme outliers (>2μm): {len(df[df['distance_to_skeleton_um'] > 2.0]):,} spots
   • Severely affected cells: {len(cell_stats[cell_stats['mean_distance'] > 1.5]):,}
   • Max spots per cell: {cell_stats['spot_count'].max():,}

🔍 Top Issues:
   • Worst cell: {cell_stats.index[0][1]} ({cell_stats.iloc[0]['mean_distance']:.3f}μm avg)
   • Most outliers in one cell: {cell_stats['spot_count'].max()} spots
"""
    
    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    # Save summary plot
    summary_file = os.path.join(output_dir, 'outlier_analysis_summary.png')
    plt.savefig(summary_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Summary visualization saved: {summary_file}")

def generate_recommendations(df, cell_stats, image_stats):
    """Generate actionable recommendations based on the analysis"""
    print("\n=== Recommendations & Action Items ===\n")
    
    # Severity assessment
    total_cells = df.groupby(['image', 'cell']).ngroups
    extreme_spots = len(df[df['distance_to_skeleton_um'] > 2.0])
    severe_cells = len(cell_stats[cell_stats['mean_distance'] > 1.5])
    
    print("🎯 Severity Assessment:")
    if extreme_spots == 0:
        print("   ✅ LOW SEVERITY: No extreme outliers (>2μm) detected")
    elif extreme_spots < 10:
        print("   ⚠️  MODERATE SEVERITY: Few extreme outliers detected")
    else:
        print("   🚨 HIGH SEVERITY: Many extreme outliers detected")
    
    # Specific recommendations
    print("\n📋 Recommended Actions:")
    
    # 1. Cell-level recommendations
    if severe_cells > 0:
        print(f"   1. 🔍 PRIORITY: Investigate {severe_cells} severely affected cells:")
        worst_cells = cell_stats.head(5)
        for i, (idx, row) in enumerate(worst_cells.iterrows(), 1):
            print(f"      • {idx[1]} in {idx[0]} - {row['spot_count']} spots, {row['mean_distance']:.3f}μm avg")
    
    # 2. Image-level recommendations  
    if len(image_stats) > 0:
        worst_image = image_stats.index[0]
        print(f"   2. 🖼️  FOCUS: Start with image '{worst_image}'")
        print(f"      • Has {image_stats.loc[worst_image, 'unique_cells']} affected cells")
        print(f"      • Average outlier distance: {image_stats.loc[worst_image, 'mean_distance']:.3f}μm")
    
    # 3. Technical recommendations
    print(f"   3. ⚙️  TECHNICAL CHECKS:")
    if 'y_flip' in cell_stats.columns:
        y_flip_rate = cell_stats['y_flip'].mean()
        if y_flip_rate > 0.5:
            print(f"      • HIGH Y-flip usage ({y_flip_rate*100:.1f}%) - check coordinate alignment")
        elif y_flip_rate < 0.1:
            print(f"      • LOW Y-flip usage ({y_flip_rate*100:.1f}%) - might need more Y-flips")
    
    print(f"      • Verify coordinate mapping for worst-performing images")
    print(f"      • Check skeleton extraction quality")
    print(f"      • Validate FISH-QUANT spot detection parameters")
    
    # 4. Threshold recommendations
    current_threshold = df['threshold_um'].iloc[0]
    q75 = np.percentile(df['distance_to_skeleton_um'], 75)
    if q75 < current_threshold * 1.5:
        print(f"   4. 📏 THRESHOLD: Consider increasing threshold to {current_threshold * 1.2:.1f}μm")
    elif df['distance_to_skeleton_um'].min() > current_threshold * 1.2:
        print(f"   4. 📏 THRESHOLD: Current threshold ({current_threshold}μm) might be too strict")
    
    print(f"\n💡 Next Steps:")
    print(f"   1. Review visualizations for patterns")
    print(f"   2. Check coordinate mapping files for worst cells") 
    print(f"   3. Validate skeleton extraction quality")
    print(f"   4. Consider adjusting alignment parameters")
    print(f"   5. Re-run analysis after corrections")

def main():
    """Main analysis function"""
    
    # File path - update this to match your file location
    file_path = "/Volumes/ExFAT/mRNA_molecules/Y333 ATP6 ATP2/interactive_batch_results/atp6_corrected_outlier_analysis_0.6_extracted_cells_conn_nonadaptive.csv"
    output_dir = "Y333 ATP6 ATP2/interactive_batch_results"
    
    print("🔬 SINGLE OUTLIER FILE ANALYSIS")
    print("=" * 50)
    
    # Load data
    df = load_outlier_data(file_path)
    if df is None:
        print("❌ Cannot proceed without valid data")
        return
    
    # Run analyses
    analyze_basic_statistics(df)
    cell_stats = analyze_by_cell(df)
    image_stats = analyze_by_image(df)
    
    # Create visualizations
    create_visualizations(df, cell_stats, image_stats, output_dir)
    
    # Generate recommendations
    generate_recommendations(df, cell_stats, image_stats)
    
    print(f"\n✅ Analysis complete! Check the output directory: {output_dir}")
    print("📊 Generated files:")
    print("   • comprehensive_outlier_analysis.png - Detailed analysis plots")
    print("   • outlier_analysis_summary.png - Quick overview")

if __name__ == "__main__":
    main() 