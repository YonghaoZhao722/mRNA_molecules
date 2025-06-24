import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import json

def parse_spots_file(spots_file):
    """Parse spots file, extract spot coordinates for each cell"""
    with open(spots_file, 'r') as f:
        lines = f.readlines()
    
    cells_data = {}
    current_cell = None
    reading_spots = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('CELL'):
            current_cell = line.split('\t')[1]
            cells_data[current_cell] = {'spots': []}
            reading_spots = False
        elif line == 'SPOTS':
            reading_spots = True
            continue
        elif reading_spots and line and not line.startswith('Pos_Y'):
            if line.startswith('CELL'):
                current_cell = line.split('\t')[1]
                cells_data[current_cell] = {'spots': []}
                reading_spots = False
            else:
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        # Coordinate format: Pos_Y, Pos_X, Pos_Z (nanometer coordinates)
                        y = float(parts[0])
                        x = float(parts[1]) 
                        z = float(parts[2])
                        cells_data[current_cell]['spots'].append([z, y, x])  # Convert to zyx format
                    except ValueError:
                        continue
    
    return cells_data

def read_skeleton_data(skeleton_file):
    """Read skeleton data"""
    df = pd.read_csv(skeleton_file, sep='\t')
    # Extract xyz coordinates
    coords = df[['x', 'y', 'z']].values
    return coords

def convert_and_scale_coordinates(spots_nm, skeleton_coords, voxel_size=(300, 160, 160)):
    """
    Convert spots coordinates and scale to skeleton coordinate system with Y-axis flip
    spots_nm: spot coordinates (z, y, x) in nanometers
    skeleton_coords: skeleton coordinates (x, y, z) in micrometers
    voxel_size: (z, y, x) in nanometers per voxel
    """
    print(f"  Original spots coordinate range (z,y,x) [nm]:")
    print(f"    Z: {spots_nm[:,0].min():.0f} - {spots_nm[:,0].max():.0f}")
    print(f"    Y: {spots_nm[:,1].min():.0f} - {spots_nm[:,1].max():.0f}") 
    print(f"    X: {spots_nm[:,2].min():.0f} - {spots_nm[:,2].max():.0f}")
    
    # Convert from nanometers to micrometers
    spots_um = spots_nm / 1000.0
    print(f"  After converting to micrometers (z,y,x) [μm]:")
    print(f"    Z: {spots_um[:,0].min():.3f} - {spots_um[:,0].max():.3f}")
    print(f"    Y: {spots_um[:,1].min():.3f} - {spots_um[:,1].max():.3f}")
    print(f"    X: {spots_um[:,2].min():.3f} - {spots_um[:,2].max():.3f}")
    
    # Rearrange to (x, y, z) format to match skeleton
    spots_um_xyz = spots_um[:, [2, 1, 0]]  # From (z, y, x) to (x, y, z)
    print(f"  After rearranging to xyz format [μm]:")
    print(f"    X: {spots_um_xyz[:,0].min():.3f} - {spots_um_xyz[:,0].max():.3f}")
    print(f"    Y: {spots_um_xyz[:,1].min():.3f} - {spots_um_xyz[:,1].max():.3f}")
    print(f"    Z: {spots_um_xyz[:,2].min():.3f} - {spots_um_xyz[:,2].max():.3f}")
    
    print(f"  Skeleton coordinate range (x,y,z) [μm]:")
    print(f"    X: {skeleton_coords[:,0].min():.3f} - {skeleton_coords[:,0].max():.3f}")
    print(f"    Y: {skeleton_coords[:,1].min():.3f} - {skeleton_coords[:,1].max():.3f}")
    print(f"    Z: {skeleton_coords[:,2].min():.3f} - {skeleton_coords[:,2].max():.3f}")
    
    # Calculate scaling factors
    if len(spots_um_xyz) > 1:
        spots_x_range = spots_um_xyz[:, 0].max() - spots_um_xyz[:, 0].min()
        spots_y_range = spots_um_xyz[:, 1].max() - spots_um_xyz[:, 1].min()
        
        skeleton_x_range = skeleton_coords[:, 0].max() - skeleton_coords[:, 0].min()
        skeleton_y_range = skeleton_coords[:, 1].max() - skeleton_coords[:, 1].min()
        
        if spots_x_range > 0 and spots_y_range > 0:
            scale_x = skeleton_x_range / spots_x_range
            scale_y = skeleton_y_range / spots_y_range
        else:
            scale_x = scale_y = 1.0
    else:
        scale_x = scale_y = 1.0
    
    print(f"  Calculated scaling factors: X={scale_x:.6f}, Y={scale_y:.6f}")
    
    # Apply Y-axis flip and scaling (default method)
    spots_flipped_y = spots_um_xyz.copy()
    spots_flipped_y[:, 1] = -spots_flipped_y[:, 1]  # Y coordinate flip
    
    spots_scaled = spots_flipped_y.copy()
    if len(spots_um_xyz) > 1:
        spots_scaled[:, 0] = (spots_scaled[:, 0] - spots_scaled[:, 0].min()) * scale_x + skeleton_coords[:, 0].min()
        spots_scaled[:, 1] = (spots_scaled[:, 1] - spots_scaled[:, 1].min()) * scale_y + skeleton_coords[:, 1].min()
    else:
        spots_scaled[:, 0] = (skeleton_coords[:, 0].min() + skeleton_coords[:, 0].max()) / 2
        spots_scaled[:, 1] = (skeleton_coords[:, 1].min() + skeleton_coords[:, 1].max()) / 2
    spots_scaled[:, 2] = spots_flipped_y[:, 2]
    
    print(f"  Y-axis flipped and scaled spots coordinate range [μm]:")
    print(f"    X: {spots_scaled[:,0].min():.3f} - {spots_scaled[:,0].max():.3f}")
    print(f"    Y: {spots_scaled[:,1].min():.3f} - {spots_scaled[:,1].max():.3f}")
    print(f"    Z: {spots_scaled[:,2].min():.3f} - {spots_scaled[:,2].max():.3f}")
    
    return {
        'aligned_spots': spots_scaled,
        'scale_x': scale_x,
        'scale_y': scale_y
    }

def calculate_min_distances(spots_coords, skeleton_coords):
    """Calculate minimum distance from each spot to skeleton"""
    if len(spots_coords) == 0 or len(skeleton_coords) == 0:
        return np.array([])
    
    # Calculate distance matrix
    distances = cdist(spots_coords, skeleton_coords)
    # Minimum distance from each spot to skeleton
    min_distances = np.min(distances, axis=1)
    
    return min_distances

def visualize_results(cell_id, spots_coords, skeleton_coords, distances, output_dir):
    """Visualize results"""
    fig = plt.figure(figsize=(15, 5))
    
    # 3D scatter plot
    ax1 = fig.add_subplot(131, projection='3d')
    if len(spots_coords) > 0:
        scatter = ax1.scatter(spots_coords[:, 0], spots_coords[:, 1], spots_coords[:, 2], 
                            c=distances, cmap='viridis', s=50, alpha=0.7)
        plt.colorbar(scatter, ax=ax1, label='Distance to skeleton (μm)')
    
    if len(skeleton_coords) > 0:
        ax1.scatter(skeleton_coords[:, 0], skeleton_coords[:, 1], skeleton_coords[:, 2], 
                   c='red', s=10, alpha=0.5, label='Skeleton')
    
    ax1.set_xlabel('X (μm)')
    ax1.set_ylabel('Y (μm)')
    ax1.set_zlabel('Z (μm)')
    ax1.set_title(f'{cell_id} - 3D View (Y-flipped)')
    ax1.legend()
    
    # XY projection
    ax2 = fig.add_subplot(132)
    if len(spots_coords) > 0:
        scatter2 = ax2.scatter(spots_coords[:, 0], spots_coords[:, 1], 
                             c=distances, cmap='viridis', s=50, alpha=0.7)
        plt.colorbar(scatter2, ax=ax2, label='Distance (μm)')
    
    if len(skeleton_coords) > 0:
        ax2.scatter(skeleton_coords[:, 0], skeleton_coords[:, 1], 
                   c='red', s=10, alpha=0.5, label='Skeleton')
    
    ax2.set_xlabel('X (μm)')
    ax2.set_ylabel('Y (μm)')
    ax2.set_title(f'{cell_id} - XY Projection (Y-flipped)')
    ax2.legend()
    ax2.axis('equal')  # Maintain aspect ratio
    
    # Distance distribution histogram
    ax3 = fig.add_subplot(133)
    if len(distances) > 0:
        ax3.hist(distances, bins=20, alpha=0.7, edgecolor='black')
        ax3.axvline(np.mean(distances), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(distances):.2f} μm')
        ax3.axvline(np.median(distances), color='orange', linestyle='--', 
                   label=f'Median: {np.median(distances):.2f} μm')
    
    ax3.set_xlabel('Distance to skeleton (μm)')
    ax3.set_ylabel('Frequency')
    ax3.set_title(f'{cell_id} - Distance Distribution')
    ax3.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{cell_id}_analysis_final.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Path settings
    spots_file = r'F:\atp\sample_atp2\yWL333_cy3_ATP2_cy5_ATP6MS2_3_w2CY3-100-_s1__spots.txt'
    mitograph_dir = r'F:\atp\sample_atp2\mitograph'
    output_dir = r'F:\atp\sample_atp2\spot_skeleton_analysis_final'
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse spots file
    print("Parsing spots file...")
    cells_data = parse_spots_file(spots_file)
    print(f"Found {len(cells_data)} cells")
    
    # Process each cell
    results_summary = []
    
    for cell_name, cell_data in cells_data.items():
        print(f"\nProcessing {cell_name}...")
        
        # Extract cell ID (assuming format is Cell_X)
        try:
            cell_id = int(cell_name.split('_')[1])
        except:
            print(f"Cannot parse cell ID: {cell_name}")
            continue
        
        # Check if corresponding skeleton file exists
        skeleton_file = os.path.join(mitograph_dir, f'cell_{cell_id:03d}.txt')
        
        if not os.path.exists(skeleton_file):
            print(f"Skeleton file not found: {skeleton_file}")
            continue
        
        # Read skeleton data
        skeleton_coords = read_skeleton_data(skeleton_file)
        
        # Get spots coordinates
        spots_nm = np.array(cell_data['spots'])
        
        if len(spots_nm) == 0:
            print(f"Cell {cell_name} has no spots")
            continue
        
        # Coordinate conversion and scaling (with Y-axis flip)
        conversion_results = convert_and_scale_coordinates(spots_nm, skeleton_coords)
        spots_aligned = conversion_results['aligned_spots']
        
        # Calculate distances
        distances = calculate_min_distances(spots_aligned, skeleton_coords)
        
        # Statistics
        stats = {
            'cell_id': cell_name,
            'transformation': 'Y-axis flip',
            'num_spots': len(spots_nm),
            'num_skeleton_points': len(skeleton_coords),
            'mean_distance': np.mean(distances),
            'median_distance': np.median(distances),
            'std_distance': np.std(distances),
            'min_distance': np.min(distances),
            'max_distance': np.max(distances),
            'scale_x': conversion_results['scale_x'],
            'scale_y': conversion_results['scale_y']
        }
        results_summary.append(stats)
        
        print(f"  Number of spots: {stats['num_spots']}")
        print(f"  Number of skeleton points: {stats['num_skeleton_points']}")
        print(f"  Mean distance: {stats['mean_distance']:.2f} μm")
        print(f"  Median distance: {stats['median_distance']:.2f} μm")
        
        # Visualization
        visualize_results(cell_name, spots_aligned, skeleton_coords, 
                         distances, output_dir)
        
        # Save detailed data
        detail_data = pd.DataFrame({
            'spot_x_um': spots_aligned[:, 0],
            'spot_y_um': spots_aligned[:, 1],
            'spot_z_um': spots_aligned[:, 2],
            'distance_to_skeleton_um': distances
        })
        detail_data.to_csv(os.path.join(output_dir, f'{cell_name}_distances_final.csv'), index=False)
    
    # Save summary results
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        summary_df.to_csv(os.path.join(output_dir, 'summary_statistics_final.csv'), index=False)
        
        # Plot summary graphs
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        axes[0, 0].bar(summary_df['cell_id'], summary_df['num_spots'])
        axes[0, 0].set_title('Number of Spots per Cell')
        axes[0, 0].set_ylabel('Number of Spots')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        axes[0, 1].bar(summary_df['cell_id'], summary_df['mean_distance'])
        axes[0, 1].set_title('Mean Distance to Skeleton (Y-flipped)')
        axes[0, 1].set_ylabel('Distance (μm)')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        axes[1, 0].bar(summary_df['cell_id'], summary_df['median_distance'])
        axes[1, 0].set_title('Median Distance to Skeleton (Y-flipped)')
        axes[1, 0].set_ylabel('Distance (μm)')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        axes[1, 1].hist(summary_df['mean_distance'], bins=10, alpha=0.7, edgecolor='black')
        axes[1, 1].set_title('Distribution of Mean Distances (Y-flipped)')
        axes[1, 1].set_xlabel('Mean Distance (μm)')
        axes[1, 1].set_ylabel('Frequency')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'summary_analysis_final.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\nAnalysis completed! Results saved in: {output_dir}")
        print(f"Processed {len(results_summary)} cells")
        
        # Print summary statistics
        print(f"\nSummary statistics:")
        print(f"  Mean distance range: {summary_df['mean_distance'].min():.2f} - {summary_df['mean_distance'].max():.2f} μm")
        print(f"  Median distance range: {summary_df['median_distance'].min():.2f} - {summary_df['median_distance'].max():.2f} μm")
        print(f"  Overall mean distance: {summary_df['mean_distance'].mean():.2f} ± {summary_df['mean_distance'].std():.2f} μm")
        
    else:
        print("No cells were successfully processed")

if __name__ == "__main__":
    main() 