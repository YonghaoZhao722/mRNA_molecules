#!/usr/bin/env python3
"""
Batch analysis of skeleton-spot alignment across all cells and image sequences
Extended from analyze_alignment.py to process the entire dataset
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import json
from pathlib import Path
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Import functions from the original analyze_alignment.py
from analyze_alignment import (
    load_cell_data_static_mapping,
    load_cell_data_dynamic_mapping,
    test_z_coordinate_corrections,
    calculate_alignment_metrics,
    analyze_z_coordinate_mapping
)

def find_all_cell_files():
    """Find all skeleton files and their corresponding spots files"""
    base_dir = "Y333 ATP6 ATP2"
    extracted_cells_dir = os.path.join(base_dir, "extracted_cells")
    spots_dir = os.path.join(base_dir, "atp6_spots")
    
    cell_files = []
    
    # Get all skeleton directories
    skeleton_dirs = glob.glob(os.path.join(extracted_cells_dir, "yWL333_*"))
    
    for skeleton_dir in skeleton_dirs:
        # Extract experiment and position info from directory name
        dir_name = os.path.basename(skeleton_dir)
        # Example: yWL333_cy3_ATP2_cy5_ATP6MS2_1_s1
        parts = dir_name.split('_')
        if len(parts) >= 6:
            experiment = parts[5]  # 1, 2, or 3
            position = parts[6]    # s1, s2, etc.
            
            # Find corresponding spots file
            spots_pattern = f"yWL333_cy3_ATP2_cy5_ATP6MS2_{experiment}_DIC_{position}_spots_*.txt"
            spots_files = glob.glob(os.path.join(spots_dir, spots_pattern))
            
            if spots_files:
                spots_file = spots_files[0]  # Take the first match
                
                # Find all cell skeleton files in this directory
                cell_skeleton_files = glob.glob(os.path.join(skeleton_dir, "cell_*.txt"))
                
                for skeleton_file in cell_skeleton_files:
                    # Extract cell number from filename
                    cell_basename = os.path.basename(skeleton_file)
                    cell_id = cell_basename.replace('cell_', '').replace('.txt', '')
                    
                    # Check if coordinate mapping file exists
                    mapping_file = os.path.join(skeleton_dir, 'coordinate_mapping.json')
                    has_mapping = os.path.exists(mapping_file)
                    
                    cell_info = {
                        'experiment': experiment,
                        'position': position,
                        'cell_id': cell_id,
                        'skeleton_file': skeleton_file,
                        'spots_file': spots_file,
                        'mapping_file': mapping_file if has_mapping else None,
                        'sequence_dir': skeleton_dir
                    }
                    cell_files.append(cell_info)
    
    print(f"Found {len(cell_files)} cells across {len(set([c['experiment'] + '_' + c['position'] for c in cell_files]))} image sequences")
    return cell_files

def analyze_single_cell(cell_info, method='static', quiet=True):
    """Analyze alignment for a single cell"""    
    try:
        skeleton_file = cell_info['skeleton_file']
        spots_file = cell_info['spots_file']
        cell_number = int(cell_info['cell_id'])
        
        # Load mapping data if available
        mapping_data = None
        if cell_info['mapping_file'] and os.path.exists(cell_info['mapping_file']):
            with open(cell_info['mapping_file'], 'r') as f:
                mapping_data = json.load(f)
        
        # Simple output suppression using devnull redirection
        if quiet:
            import sys
            old_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
        
        try:
            # Load cell data using specified method
            if method == 'static':
                spots_coords, skeleton_coords, pixel_xy, pixel_z = load_cell_data_static_mapping(
                    skeleton_file, spots_file, cell_number, mapping_data
                )
            else:  # dynamic
                spots_coords, skeleton_coords, pixel_xy, pixel_z = load_cell_data_dynamic_mapping(
                    skeleton_file, spots_file, cell_number, mapping_data
                )
            
            # Calculate basic alignment metrics
            metrics = calculate_alignment_metrics(spots_coords, skeleton_coords, method_name="")
            
            # Test Z-coordinate corrections
            spots_z_corrected, optimal_z_offset, z_correction_results = test_z_coordinate_corrections(
                skeleton_coords, spots_coords
            )
            
            # Calculate metrics with Z-correction
            corrected_metrics = calculate_alignment_metrics(spots_z_corrected, skeleton_coords, method_name="")
            
            # Analyze Z-coordinate mapping
            z_differences = analyze_z_coordinate_mapping(skeleton_coords, spots_coords)
            
        finally:
            # Restore stdout
            if quiet:
                sys.stdout.close()
                sys.stdout = old_stdout
        
        results = {
            'experiment': cell_info['experiment'],
            'position': cell_info['position'], 
            'cell_id': cell_info['cell_id'],
            'method': method,
            'n_spots': len(spots_coords),
            'n_skeleton_points': len(skeleton_coords),
            'has_mapping': cell_info['mapping_file'] is not None,
            'pixel_xy': pixel_xy,
            'pixel_z': pixel_z,
            # Original metrics
            'center_distance': metrics['center_distance'],
            'mean_distance': metrics['mean_distance'],
            'median_distance': metrics['median_distance'],
            'std_distance': metrics['std_distance'],
            'within_1um': metrics['within_1um'],
            'within_0_5um': metrics['within_0_5um'],
            'within_0_1um': metrics['within_0_1um'],
            # Z-correction metrics
            'optimal_z_offset': optimal_z_offset,
            'corrected_mean_distance': corrected_metrics['mean_distance'],
            'corrected_within_0_5um': corrected_metrics['within_0_5um'],
            'z_diff_mean': np.mean(z_differences),
            'z_diff_std': np.std(z_differences),
            'z_diff_median': np.median(z_differences),
            'skeleton_file': skeleton_file,
            'spots_file': spots_file
        }
        
        return results
        
    except Exception as e:
        print(f"Error processing {cell_info['experiment']}_{cell_info['position']}_cell_{cell_info['cell_id']}: {e}")
        return None

def batch_analyze_alignment(method='static', max_cells=None):
    """Perform batch alignment analysis on all cells"""
    print(f"=== BATCH ALIGNMENT ANALYSIS ({method.upper()} METHOD) ===")
    
    # Find all cell files
    cell_files = find_all_cell_files()
    
    if max_cells:
        cell_files = cell_files[:max_cells]
        print(f"Limiting analysis to first {max_cells} cells")
    
    # Analyze each cell
    all_results = []
    failed_cells = []
    
    print(f"Analyzing {len(cell_files)} cells...")
    
    for cell_info in tqdm(cell_files, desc="Processing cells"):
        result = analyze_single_cell(cell_info, method=method, quiet=True)
        if result:
            all_results.append(result)
        else:
            failed_cells.append(f"{cell_info['experiment']}_{cell_info['position']}_cell_{cell_info['cell_id']}")
    
    print(f"Successfully analyzed {len(all_results)} cells")
    if failed_cells:
        print(f"Failed to analyze {len(failed_cells)} cells: {failed_cells[:10]}{'...' if len(failed_cells) > 10 else ''}")
    
    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    
    return df, failed_cells

def create_comprehensive_plots(df, method='static'):
    """Create comprehensive visualization of batch analysis results"""
    
    # Check if DataFrame is empty
    if len(df) == 0:
        print(f"Warning: No data to plot for {method} method")
        return None
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(3, 4, figsize=(24, 18))
    fig.suptitle(f'Comprehensive Alignment Analysis - {method.upper()} Method\n'
                f'Total Cells: {len(df)}', fontsize=16, fontweight='bold')
    
    # 1. Mean distance distribution
    axes[0,0].hist(df['mean_distance'], bins=50, alpha=0.7, edgecolor='black')
    axes[0,0].set_xlabel('Mean Distance to Skeleton (μm)')
    axes[0,0].set_ylabel('Number of Cells')
    axes[0,0].set_title('Distribution of Mean Distances')
    axes[0,0].axvline(df['mean_distance'].mean(), color='red', linestyle='--', 
                     label=f'Overall Mean: {df["mean_distance"].mean():.3f}μm')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # 2. Z-offset distribution  
    axes[0,1].hist(df['optimal_z_offset'], bins=50, alpha=0.7, edgecolor='black', color='orange')
    axes[0,1].set_xlabel('Optimal Z-Offset (μm)')
    axes[0,1].set_ylabel('Number of Cells')
    axes[0,1].set_title('Distribution of Optimal Z-Offsets')
    axes[0,1].axvline(df['optimal_z_offset'].mean(), color='red', linestyle='--',
                     label=f'Mean: {df["optimal_z_offset"].mean():.3f}μm')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # 3. Improvement from Z-correction
    improvement = df['mean_distance'] - df['corrected_mean_distance']
    axes[0,2].hist(improvement, bins=50, alpha=0.7, edgecolor='black', color='green')
    axes[0,2].set_xlabel('Distance Improvement (μm)')
    axes[0,2].set_ylabel('Number of Cells')
    axes[0,2].set_title('Improvement from Z-Correction')
    axes[0,2].axvline(improvement.mean(), color='red', linestyle='--',
                     label=f'Mean Improvement: {improvement.mean():.3f}μm')
    axes[0,2].legend()
    axes[0,2].grid(True, alpha=0.3)
    
    # 4. Within 0.5μm percentage
    axes[0,3].hist(df['within_0_5um'] * 100, bins=30, alpha=0.7, edgecolor='black', color='purple')
    axes[0,3].set_xlabel('Spots Within 0.5μm (%)')
    axes[0,3].set_ylabel('Number of Cells')
    axes[0,3].set_title('Percentage of Well-Aligned Spots')
    axes[0,3].axvline((df['within_0_5um'] * 100).mean(), color='red', linestyle='--',
                     label=f'Mean: {(df["within_0_5um"] * 100).mean():.1f}%')
    axes[0,3].legend()
    axes[0,3].grid(True, alpha=0.3)
    
    # 5. By experiment
    exp_data = df.groupby('experiment').agg({
        'mean_distance': 'mean',
        'optimal_z_offset': 'mean',
        'within_0_5um': 'mean'
    }).reset_index()
    
    x_pos = np.arange(len(exp_data))
    axes[1,0].bar(x_pos, exp_data['mean_distance'], alpha=0.7)
    axes[1,0].set_xlabel('Experiment')
    axes[1,0].set_ylabel('Mean Distance (μm)')
    axes[1,0].set_title('Mean Distance by Experiment')
    axes[1,0].set_xticks(x_pos)
    axes[1,0].set_xticklabels(exp_data['experiment'])
    axes[1,0].grid(True, alpha=0.3)
    
    # 6. Z-offset by experiment
    axes[1,1].bar(x_pos, exp_data['optimal_z_offset'], alpha=0.7, color='orange')
    axes[1,1].set_xlabel('Experiment')
    axes[1,1].set_ylabel('Mean Z-Offset (μm)')
    axes[1,1].set_title('Z-Offset by Experiment')
    axes[1,1].set_xticks(x_pos)
    axes[1,1].set_xticklabels(exp_data['experiment'])
    axes[1,1].grid(True, alpha=0.3)
    
    # 7. Correlation: Z-offset vs improvement
    axes[1,2].scatter(df['optimal_z_offset'], improvement, alpha=0.6)
    axes[1,2].set_xlabel('Optimal Z-Offset (μm)')
    axes[1,2].set_ylabel('Distance Improvement (μm)')
    axes[1,2].set_title('Z-Offset vs Improvement')
    axes[1,2].grid(True, alpha=0.3)
    
    # Add correlation coefficient
    corr = np.corrcoef(df['optimal_z_offset'], improvement)[0,1]
    axes[1,2].text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=axes[1,2].transAxes,
                  bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 8. Number of spots vs alignment quality
    axes[1,3].scatter(df['n_spots'], df['mean_distance'], alpha=0.6, color='purple')
    axes[1,3].set_xlabel('Number of Spots')
    axes[1,3].set_ylabel('Mean Distance (μm)')
    axes[1,3].set_title('Spot Count vs Alignment Quality')
    axes[1,3].grid(True, alpha=0.3)
    
    # 9. Position effects (if multiple positions)
    if 'position' in df.columns and len(df['position'].unique()) > 1:
        pos_data = df.groupby('position')['mean_distance'].mean().reset_index()
        pos_data['pos_num'] = pos_data['position'].str.extract('(\d+)').astype(int)
        pos_data = pos_data.sort_values('pos_num')
        
        axes[2,0].plot(pos_data['pos_num'], pos_data['mean_distance'], 'o-', alpha=0.7)
        axes[2,0].set_xlabel('Position Number')
        axes[2,0].set_ylabel('Mean Distance (μm)')
        axes[2,0].set_title('Alignment Quality by Position')
        axes[2,0].grid(True, alpha=0.3)
    else:
        axes[2,0].text(0.5, 0.5, 'Position analysis\nnot applicable', 
                      ha='center', va='center', transform=axes[2,0].transAxes)
    
    # 10. Z-difference analysis
    axes[2,1].hist(df['z_diff_mean'], bins=50, alpha=0.7, edgecolor='black', color='brown')
    axes[2,1].set_xlabel('Mean Z-Difference (μm)')
    axes[2,1].set_ylabel('Number of Cells')
    axes[2,1].set_title('Systematic Z-Difference Distribution')
    axes[2,1].axvline(df['z_diff_mean'].mean(), color='red', linestyle='--',
                     label=f'Overall Mean: {df["z_diff_mean"].mean():.3f}μm')
    axes[2,1].legend()
    axes[2,1].grid(True, alpha=0.3)
    
    # 11. Before vs After Z-correction scatter
    axes[2,2].scatter(df['mean_distance'], df['corrected_mean_distance'], alpha=0.6, color='teal')
    axes[2,2].plot([0, df['mean_distance'].max()], [0, df['mean_distance'].max()], 'r--', alpha=0.8)
    axes[2,2].set_xlabel('Original Mean Distance (μm)')
    axes[2,2].set_ylabel('Z-Corrected Mean Distance (μm)')
    axes[2,2].set_title('Before vs After Z-Correction')
    axes[2,2].grid(True, alpha=0.3)
    
    # 12. Summary statistics text
    axes[2,3].axis('off')
    summary_stats = f"""BATCH ANALYSIS SUMMARY
    
Total Cells Analyzed: {len(df)}
Experiments: {len(df['experiment'].unique())}
Positions: {len(df['position'].unique())}

ALIGNMENT METRICS:
Mean Distance: {df['mean_distance'].mean():.3f} ± {df['mean_distance'].std():.3f} μm
Median Distance: {df['median_distance'].mean():.3f} μm
Within 0.5μm: {(df['within_0_5um'] * 100).mean():.1f}%

Z-CORRECTION ANALYSIS:
Mean Z-Offset: {df['optimal_z_offset'].mean():.3f} ± {df['optimal_z_offset'].std():.3f} μm
Mean Improvement: {improvement.mean():.3f} μm
Cells Needing Correction (>50nm): {sum(abs(df['optimal_z_offset']) > 0.05)} ({sum(abs(df['optimal_z_offset']) > 0.05)/len(df)*100:.1f}%)

SYSTEMATIC ISSUES:
Mean Z-Difference: {df['z_diff_mean'].mean():.3f} μm
Z-Offset Range: [{df['optimal_z_offset'].min():.3f}, {df['optimal_z_offset'].max():.3f}] μm
"""
    
    axes[2,3].text(0.05, 0.95, summary_stats, transform=axes[2,3].transAxes,
                  fontsize=10, verticalalignment='top', fontfamily='monospace',
                  bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    output_filename = f'batch_alignment_analysis_{method}.png'
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Comprehensive analysis plot saved to: {output_filename}")
    
    return fig

def save_results(df, method='static'):
    """Save detailed results to CSV files"""
    
    # Save full results
    output_csv = f'batch_alignment_results_{method}.csv'
    df.to_csv(output_csv, index=False)
    print(f"Detailed results saved to: {output_csv}")
    
    # Check if DataFrame is empty
    if len(df) == 0:
        print(f"Warning: No successful analyses to summarize for {method} method")
        summary_stats = {
            'total_cells': 0,
            'mean_distance_mean': None,
            'mean_distance_std': None,
            'optimal_z_offset_mean': None,
            'optimal_z_offset_std': None,
            'within_0_5um_mean': None,
            'improvement_mean': None,
            'cells_needing_z_correction': 0,
            'systematic_z_bias': None
        }
    else:
        # Save summary statistics
        summary_stats = {
            'total_cells': len(df),
            'mean_distance_mean': df['mean_distance'].mean(),
            'mean_distance_std': df['mean_distance'].std(),
            'optimal_z_offset_mean': df['optimal_z_offset'].mean(),
            'optimal_z_offset_std': df['optimal_z_offset'].std(),
            'within_0_5um_mean': df['within_0_5um'].mean(),
            'improvement_mean': (df['mean_distance'] - df['corrected_mean_distance']).mean(),
            'cells_needing_z_correction': sum(abs(df['optimal_z_offset']) > 0.05),
            'systematic_z_bias': df['z_diff_mean'].mean()
        }
    
    summary_df = pd.DataFrame([summary_stats])
    summary_csv = f'batch_alignment_summary_{method}.csv'
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary statistics saved to: {summary_csv}")
    
    return output_csv, summary_csv

def compare_methods():
    """Compare static vs dynamic methods on the same dataset"""
    print("=== COMPARING STATIC vs DYNAMIC METHODS ===")
    
    # Analyze with both methods (limit to subset for comparison)
    print("Analyzing with STATIC method...")
    df_static, _ = batch_analyze_alignment(method='static', max_cells=100)
    
    print("Analyzing with DYNAMIC method...")
    df_dynamic, _ = batch_analyze_alignment(method='dynamic', max_cells=100)
    
    # Merge results for comparison
    df_static['method'] = 'Static'
    df_dynamic['method'] = 'Dynamic'
    df_combined = pd.concat([df_static, df_dynamic], ignore_index=True)
    
    # Create comparison plot
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Static vs Dynamic Method Comparison', fontsize=16, fontweight='bold')
    
    # Mean distance comparison
    sns.boxplot(data=df_combined, x='method', y='mean_distance', ax=axes[0,0])
    axes[0,0].set_title('Mean Distance Distribution')
    axes[0,0].set_ylabel('Mean Distance (μm)')
    
    # Within 0.5μm comparison
    sns.boxplot(data=df_combined, x='method', y='within_0_5um', ax=axes[0,1])
    axes[0,1].set_title('Within 0.5μm Percentage')
    axes[0,1].set_ylabel('Percentage Within 0.5μm')
    
    # Z-offset comparison
    sns.boxplot(data=df_combined, x='method', y='optimal_z_offset', ax=axes[0,2])
    axes[0,2].set_title('Optimal Z-Offset Distribution')
    axes[0,2].set_ylabel('Z-Offset (μm)')
    
    # Experiment-wise comparison
    exp_comparison = df_combined.groupby(['experiment', 'method'])['mean_distance'].mean().reset_index()
    sns.barplot(data=exp_comparison, x='experiment', y='mean_distance', hue='method', ax=axes[1,0])
    axes[1,0].set_title('Mean Distance by Experiment')
    axes[1,0].set_ylabel('Mean Distance (μm)')
    
    # Statistical comparison text
    axes[1,1].axis('off')
    from scipy import stats
    static_distances = df_static['mean_distance']
    dynamic_distances = df_dynamic['mean_distance']
    t_stat, p_value = stats.ttest_ind(static_distances, dynamic_distances)
    
    comparison_text = f"""STATISTICAL COMPARISON

Static Method:
  Mean Distance: {static_distances.mean():.3f} ± {static_distances.std():.3f} μm
  Within 0.5μm: {(df_static['within_0_5um'].mean() * 100):.1f}%

Dynamic Method:
  Mean Distance: {dynamic_distances.mean():.3f} ± {dynamic_distances.std():.3f} μm
  Within 0.5μm: {(df_dynamic['within_0_5um'].mean() * 100):.1f}%

T-test Results:
  t-statistic: {t_stat:.3f}
  p-value: {p_value:.3e}
  {'Significant difference' if p_value < 0.05 else 'No significant difference'}

RECOMMENDATION:
{'Static method preferred' if static_distances.mean() < dynamic_distances.mean() else 'Dynamic method preferred'}
based on mean distance metric.
"""
    
    axes[1,1].text(0.05, 0.95, comparison_text, transform=axes[1,1].transAxes,
                  fontsize=11, verticalalignment='top', fontfamily='monospace',
                  bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # Hide last subplot
    axes[1,2].axis('off')
    
    plt.tight_layout()
    plt.savefig('static_vs_dynamic_comparison.png', dpi=150, bbox_inches='tight')
    print("Method comparison plot saved to: static_vs_dynamic_comparison.png")
    
    return df_combined

def main():
    """Main function to run batch analysis"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Batch alignment analysis for all cells')
    parser.add_argument('--method', choices=['static', 'dynamic', 'both'], default='static',
                       help='Analysis method to use')
    parser.add_argument('--max_cells', type=int, default=None,
                       help='Maximum number of cells to analyze (for testing)')
    parser.add_argument('--compare', action='store_true',
                       help='Compare static vs dynamic methods')
    
    args = parser.parse_args()
    
    if args.compare:
        df_combined = compare_methods()
        return
    
    if args.method == 'both':
        # Analyze with both methods
        print("Analyzing with STATIC method...")
        df_static, failed_static = batch_analyze_alignment(method='static', max_cells=args.max_cells)
        save_results(df_static, method='static')
        if len(df_static) > 0:
            create_comprehensive_plots(df_static, method='static')
        
        print("\nAnalyzing with DYNAMIC method...")
        df_dynamic, failed_dynamic = batch_analyze_alignment(method='dynamic', max_cells=args.max_cells)
        save_results(df_dynamic, method='dynamic') 
        if len(df_dynamic) > 0:
            create_comprehensive_plots(df_dynamic, method='dynamic')
        
    else:
        # Analyze with single method
        df, failed = batch_analyze_alignment(method=args.method, max_cells=args.max_cells)
        save_results(df, method=args.method)
        if len(df) > 0:
            create_comprehensive_plots(df, method=args.method)
        else:
            print(f"Skipping plot generation due to no successful analyses")
    
    print("\n=== BATCH ANALYSIS COMPLETE ===")

if __name__ == "__main__":
    main() 