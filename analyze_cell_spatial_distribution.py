#!/usr/bin/env python3
"""
Analyze spatial distribution of skeleton and spots for each cell
Calculate maximum distances (ranges) in x, y, z directions and generate histograms
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import tifffile
from skimage.measure import find_contours

class CellSpatialDistributionAnalyzer:
    def __init__(self, base_dir='Y333 ATP6 ATP2'):
        self.base_dir = base_dir
        self.skeleton_type = 'extracted_cells_30'
        self.skeleton_root = os.path.join(self.base_dir, self.skeleton_type)
        self.channel = 'atp6_corrected'
        self.spots_root = os.path.join(self.base_dir, f'{self.channel}_spots')
        
        # Data storage
        self.available_images = {}
        self.coordinate_mappings = {}
        self.analysis_results = []
        
        # Scan and load data
        self.scan_available_images_and_cells()
    
    def load_spots_fishquant_analyze_method(self, file_path, cell_number=1, flip_y=True, mapping_data=None, silent=True):
        """Load FISH-QUANT spots data for specified cell"""
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find pixel size information
        pixel_xy, pixel_z = None, None
        for i, line in enumerate(lines):
            if line.startswith('Pix-XY'):
                if i + 1 < len(lines):
                    next_line = lines[i + 1].strip()
                    parts = next_line.split('\t')
                    if len(parts) >= 2:
                        pixel_xy = float(parts[0])
                        pixel_z = float(parts[1])
                break
        
        if not silent:
            print(f"FISH-QUANT pixel size: XY={pixel_xy}nm, Z={pixel_z}nm")
        
        # Find SPOTS data for the specified cell
        target_cell = f"CELL\tCell_{cell_number}"
        cell_found = False
        spots_start = -1
        spots_end = -1
        
        for i, line in enumerate(lines):
            line_stripped = line.strip()
            
            if line_stripped == target_cell:
                cell_found = True
                if not silent:
                    print(f"Found target cell: Cell_{cell_number}")
                continue
            
            if cell_found and line.startswith('Pos_Y'):
                spots_start = i
                continue
            
            if cell_found and spots_start != -1 and line.startswith('CELL'):
                spots_end = i
                break
                
            if cell_found and spots_start == -1 and line.startswith('CELL'):
                if not silent:
                    print(f"Cell_{cell_number} has no spots data")
                raise ValueError(f"Cell_{cell_number} has no spots data")
        
        if cell_found and spots_start != -1 and spots_end == -1:
            spots_end = len(lines)
        
        if not cell_found:
            raise ValueError(f"Cell_{cell_number} data not found")
        
        if spots_start == -1:
            raise ValueError(f"Cell_{cell_number} has no spots data")
        
        # Read spots data for this cell
        spots_data = []
        for i in range(spots_start + 1, spots_end):
            line = lines[i].strip()
            if line and not line.startswith('X_POS') and not line.startswith('Y_POS'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        y, x, z = float(parts[0]), float(parts[1]), float(parts[2])
                        spots_data.append([y, x, z])
                    except ValueError:
                        continue
        
        coords = np.array(spots_data)
        
        # Apply Y-axis flip if needed
        if flip_y and mapping_data:
            cell_name = f"cell_{cell_number:03d}"
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                y_start = crop_info['y_start']
                y_end = crop_info['y_end']
                
                y_start_nm = y_start * pixel_xy
                y_end_nm = y_end * pixel_xy
                
                flip_center = y_start_nm + y_end_nm
                coords[:, 0] = flip_center - coords[:, 0]
                
                if not silent:
                    print(f"Y-axis flip applied")
        
        if not silent:
            print(f"Cell_{cell_number} spots data loaded: {len(coords)} points")
        
        return coords, pixel_xy, pixel_z

    def load_skeleton_txt_analyze_method(self, file_path, mapping_data=None, cell_name="cell_001", pixel_size_xy=0.0645, silent=True):
        """Load skeleton data and convert to absolute coordinates"""
        df = pd.read_csv(file_path, sep='\t')
        coords = df[['x', 'y', 'z']].values
        
        if not silent:
            print(f"Skeleton data loaded: {len(coords)} points")
        
        # Convert to absolute coordinates if mapping data is provided
        if mapping_data and cell_name in mapping_data:
            crop_info = mapping_data[cell_name]['crop_region']
            x_offset = crop_info['x_offset']
            y_offset = crop_info['y_offset']
            
            offset_x_um = x_offset * pixel_size_xy
            offset_y_um = y_offset * pixel_size_xy
            
            coords[:, 0] += offset_x_um
            coords[:, 1] += offset_y_um
            
            if not silent:
                print(f"Applied coordinate offset: X+{offset_x_um:.3f}μm, Y+{offset_y_um:.3f}μm")
        
        return coords

    def scan_available_images_and_cells(self):
        """Scan all images and cells"""
        self.available_images = {}
        if not os.path.exists(self.skeleton_root):
            print(f"Skeleton root not found: {self.skeleton_root}")
            return
            
        spot_files = [f for f in os.listdir(self.spots_root) if f.endswith('.txt') and '_spots' in f]
        
        for image_folder in os.listdir(self.skeleton_root):
            image_path = os.path.join(self.skeleton_root, image_folder)
            if not os.path.isdir(image_path):
                continue

            # Load coordinate mapping
            mapping_file = os.path.join(image_path, 'coordinate_mapping.json')
            if os.path.exists(mapping_file):
                with open(mapping_file, 'r') as f:
                    self.coordinate_mappings[image_folder] = json.load(f)
            else:
                print(f"Warning: coordinate_mapping.json not found in {image_path}")
                continue

            # Match spot file
            parts = image_folder.split('_')
            if len(parts) < 2:
                continue
            index = parts[-2]
            fov = parts[-1]
            
            matched_spot = None
            for spot_file in spot_files:
                if f'_{index}_' in spot_file and f'_{fov}_' in spot_file and '_spots' in spot_file:
                    matched_spot = os.path.join(self.spots_root, spot_file)
                    break
            
            if not matched_spot:
                print(f"Spot file not found for {image_folder}")
                continue

            # Scan skeleton files
            skeleton_files = [f for f in os.listdir(image_path) if f.startswith('cell_') and f.endswith('.txt')]
            image_cells = {}
            
            for skel_file in skeleton_files:
                try:
                    cell_num = int(skel_file.split('_')[1].split('.')[0])
                    cell_name = f"Cell_{cell_num}"
                    
                    # Verify spots data exists
                    try:
                        test_coords, _, _ = self.load_spots_fishquant_analyze_method(
                            matched_spot, 
                            cell_number=cell_num, 
                            flip_y=False,
                            mapping_data=None,
                            silent=True
                        )
                        
                        if len(test_coords) > 0:
                            image_cells[cell_name] = {
                                'spots_file': matched_spot,
                                'skeleton_file': os.path.join(image_path, skel_file)
                            }
                            print(f"Found {image_folder}-{cell_name} ({len(test_coords)} spots)")
                            
                    except ValueError:
                        continue
                        
                except Exception as e:
                    print(f"Error parsing {skel_file}: {e}")
            
            if image_cells:
                self.available_images[image_folder] = image_cells
        
        print(f"Total images found: {len(self.available_images)}")
        total_cells = sum(len(cells) for cells in self.available_images.values())
        print(f"Total cells found: {total_cells}")

    def calculate_spatial_ranges(self, coords):
        """Calculate spatial ranges (max - min) for x, y, z coordinates"""
        if len(coords) == 0:
            return 0.0, 0.0, 0.0
        
        x_range = coords[:, 0].max() - coords[:, 0].min()
        y_range = coords[:, 1].max() - coords[:, 1].min()
        z_range = coords[:, 2].max() - coords[:, 2].min()
        
        return x_range, y_range, z_range

    def analyze_all_cells(self):
        """Analyze spatial distribution for all cells"""
        print("Starting spatial distribution analysis...")
        skipped_cells = []
        
        for image_name, cells in self.available_images.items():
            for cell_name, cell_info in cells.items():
                try:
                    print(f"Processing {image_name} - {cell_name}")
                    
                    # Get cell information
                    cell_number = int(cell_name.split("_")[1])
                    cell_id_str = f"{cell_number:03d}"
                    mapping_cell_name = f'cell_{cell_id_str}'
                    mapping_data = self.coordinate_mappings.get(image_name, {})
                    
                    # Load spots data
                    spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                        cell_info['spots_file'], 
                        cell_number=cell_number, 
                        flip_y=True, 
                        mapping_data=mapping_data,
                        silent=True
                    )
                    
                    # Load skeleton data
                    skeleton_coords = self.load_skeleton_txt_analyze_method(
                        cell_info['skeleton_file'], 
                        mapping_data=mapping_data,
                        cell_name=mapping_cell_name,
                        pixel_size_xy=pixel_xy/1000,
                        silent=True
                    )
                    
                    # Convert spots coordinates
                    spots_nm_xyz = spots_coords[:, [1, 0, 2]]  # (y,x,z) to (x,y,z)
                    spots_um_xyz = spots_nm_xyz / 1000.0  # Convert to micrometers
                    
                    # Filter out cells with no spots or no skeleton
                    if len(skeleton_coords) == 0:
                        print(f"  Skipping {image_name}-{cell_name}: No skeleton data")
                        skipped_cells.append(f"{image_name}-{cell_name}: No skeleton data")
                        continue
                    
                    if len(spots_um_xyz) == 0:
                        print(f"  Skipping {image_name}-{cell_name}: No spots data")
                        skipped_cells.append(f"{image_name}-{cell_name}: No spots data")
                        continue
                    
                    # Calculate spatial ranges
                    skeleton_x_range, skeleton_y_range, skeleton_z_range = self.calculate_spatial_ranges(skeleton_coords)
                    spots_x_range, spots_y_range, spots_z_range = self.calculate_spatial_ranges(spots_um_xyz)
                    
                    # Filter out cells with zero ranges (all points at same location)
                    skeleton_max_range = max(skeleton_x_range, skeleton_y_range, skeleton_z_range)
                    spots_max_range = max(spots_x_range, spots_y_range, spots_z_range)
                    
                    if skeleton_max_range == 0:
                        print(f"  Skipping {image_name}-{cell_name}: Skeleton has zero range (all points at same location)")
                        skipped_cells.append(f"{image_name}-{cell_name}: Skeleton zero range")
                        continue
                        
                    if spots_max_range == 0:
                        print(f"  Skipping {image_name}-{cell_name}: Spots have zero range (all points at same location)")
                        skipped_cells.append(f"{image_name}-{cell_name}: Spots zero range")
                        continue
                    
                    print(f"  Valid cell: {len(skeleton_coords)} skeleton points, {len(spots_um_xyz)} spots")
                    print(f"    Skeleton ranges: X={skeleton_x_range:.3f}, Y={skeleton_y_range:.3f}, Z={skeleton_z_range:.3f}")
                    print(f"    Spots ranges: X={spots_x_range:.3f}, Y={spots_y_range:.3f}, Z={spots_z_range:.3f}")
                    
                    # Combined data ranges
                    if len(skeleton_coords) > 0 and len(spots_um_xyz) > 0:
                        combined_coords = np.vstack([skeleton_coords, spots_um_xyz])
                    elif len(skeleton_coords) > 0:
                        combined_coords = skeleton_coords
                    else:
                        combined_coords = spots_um_xyz
                    
                    combined_x_range, combined_y_range, combined_z_range = self.calculate_spatial_ranges(combined_coords)
                    
                    # Store results
                    result = {
                        'image_name': image_name,
                        'cell_name': cell_name,
                        'num_skeleton_points': len(skeleton_coords),
                        'num_spot_points': len(spots_um_xyz),
                        'skeleton_x_range_um': skeleton_x_range,
                        'skeleton_y_range_um': skeleton_y_range,
                        'skeleton_z_range_um': skeleton_z_range,
                        'spots_x_range_um': spots_x_range,
                        'spots_y_range_um': spots_y_range,
                        'spots_z_range_um': spots_z_range,
                        'combined_x_range_um': combined_x_range,
                        'combined_y_range_um': combined_y_range,
                        'combined_z_range_um': combined_z_range,
                        'pixel_size_xy_nm': pixel_xy,
                        'pixel_size_z_nm': pixel_z
                    }
                    
                    self.analysis_results.append(result)
                    
                except Exception as e:
                    print(f"Error processing {image_name}-{cell_name}: {e}")
                    skipped_cells.append(f"{image_name}-{cell_name}: Error - {e}")
                    continue
        
        print(f"\nAnalysis completed:")
        print(f"  Valid cells analyzed: {len(self.analysis_results)}")
        print(f"  Cells skipped: {len(skipped_cells)}")
        
        if skipped_cells:
            print("\nSkipped cells summary:")
            for skipped in skipped_cells:
                print(f"  - {skipped}")

    def create_histograms(self):
        """Create histograms for all spatial range measurements"""
        if not self.analysis_results:
            print("No analysis results to plot")
            return
        
        df = pd.DataFrame(self.analysis_results)
        
        # Define the measurements to plot
        measurements = [
            ('skeleton_x_range_um', 'Skeleton X Range'),
            ('skeleton_y_range_um', 'Skeleton Y Range'), 
            ('skeleton_z_range_um', 'Skeleton Z Range'),
            ('spots_x_range_um', 'Spots X Range'),
            ('spots_y_range_um', 'Spots Y Range'),
            ('spots_z_range_um', 'Spots Z Range'),
            ('combined_x_range_um', 'Combined X Range'),
            ('combined_y_range_um', 'Combined Y Range'),
            ('combined_z_range_um', 'Combined Z Range')
        ]
        
        # Create output directory
        output_dir = os.path.join(self.base_dir, 'spatial_distribution_analysis')
        os.makedirs(output_dir, exist_ok=True)
        
        # Create histograms for each measurement
        for col_name, display_name in measurements:
            data = df[col_name].dropna()
            if len(data) == 0:
                continue
                
            # Create figure with subplots
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Regular histogram
            ax1.hist(data, bins=30, alpha=0.7, edgecolor='black', color='skyblue')
            ax1.set_xlabel('Range (μm)')
            ax1.set_ylabel('Number of Cells')
            ax1.set_title(f'{display_name} - Linear Scale')
            ax1.grid(True, alpha=0.3)
            
            # Add statistics
            stats_text = f'N={len(data)}\nMean: {data.mean():.3f}\nMedian: {data.median():.3f}\nStd: {data.std():.3f}'
            ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes, 
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            # Log+1 histogram (handles zero values)
            log_data = np.log10(data + 1)
            ax2.hist(log_data, bins=30, alpha=0.7, edgecolor='black', color='lightcoral')
            ax2.set_xlabel('Log10(Range + 1) (μm)')
            ax2.set_ylabel('Number of Cells')
            ax2.set_title(f'{display_name} - Log10(x+1) Scale')
            ax2.grid(True, alpha=0.3)
            
            # Add log statistics
            log_stats_text = f'N={len(data)}\nLog+1 Mean: {log_data.mean():.3f}\nLog+1 Std: {log_data.std():.3f}'
            ax2.text(0.98, 0.98, log_stats_text, transform=ax2.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            plt.tight_layout()
            
            # Save figure
            filename = f'{col_name}_histogram.png'
            filepath = os.path.join(output_dir, filename)
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Saved histogram: {filename}")
        
        # Save raw data
        csv_file = os.path.join(output_dir, 'spatial_distribution_data.csv')
        df.to_csv(csv_file, index=False)
        print(f"Saved data to: {csv_file}")
        
        # Create summary statistics
        self.create_summary_report(df, output_dir)

    def create_summary_report(self, df, output_dir):
        """Create a summary report of the analysis"""
        summary_file = os.path.join(output_dir, 'analysis_summary.txt')
        
        with open(summary_file, 'w') as f:
            f.write("SPATIAL DISTRIBUTION ANALYSIS SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total cells analyzed: {len(df)}\n")
            f.write(f"Channel: {self.channel}\n")
            f.write(f"Skeleton type: {self.skeleton_type}\n\n")
            
            # Summary statistics for each measurement
            measurements = [
                ('skeleton_x_range_um', 'Skeleton X Range'),
                ('skeleton_y_range_um', 'Skeleton Y Range'), 
                ('skeleton_z_range_um', 'Skeleton Z Range'),
                ('spots_x_range_um', 'Spots X Range'),
                ('spots_y_range_um', 'Spots Y Range'),
                ('spots_z_range_um', 'Spots Z Range'),
                ('combined_x_range_um', 'Combined X Range'),
                ('combined_y_range_um', 'Combined Y Range'),
                ('combined_z_range_um', 'Combined Z Range')
            ]
            
            for col_name, display_name in measurements:
                data = df[col_name].dropna()
                if len(data) > 0:
                    f.write(f"{display_name}:\n")
                    f.write(f"  Count: {len(data)}\n")
                    f.write(f"  Mean: {data.mean():.4f} μm\n")
                    f.write(f"  Median: {data.median():.4f} μm\n")
                    f.write(f"  Std: {data.std():.4f} μm\n")
                    f.write(f"  Min: {data.min():.4f} μm\n")
                    f.write(f"  Max: {data.max():.4f} μm\n")
                    f.write(f"  25th percentile: {data.quantile(0.25):.4f} μm\n")
                    f.write(f"  75th percentile: {data.quantile(0.75):.4f} μm\n\n")
            
            # Data point counts
            f.write("DATA POINT STATISTICS:\n")
            f.write(f"Average skeleton points per cell: {df['num_skeleton_points'].mean():.1f}\n")
            f.write(f"Average spot points per cell: {df['num_spot_points'].mean():.1f}\n")
            f.write(f"Total skeleton points: {df['num_skeleton_points'].sum()}\n")
            f.write(f"Total spot points: {df['num_spot_points'].sum()}\n")
        
        print(f"Saved summary report: analysis_summary.txt")

def main():
    # Create analyzer
    analyzer = CellSpatialDistributionAnalyzer()
    
    # Run analysis
    analyzer.analyze_all_cells()
    
    # Create histograms
    analyzer.create_histograms()
    
    print("Analysis completed!")

if __name__ == "__main__":
    main() 