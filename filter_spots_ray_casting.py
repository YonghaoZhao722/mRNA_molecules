#!/usr/bin/env python3
"""
Script to filter spot data using ray casting algorithm.
Removes spots that are outside the cell outline polygon.

Usage: python filter_spots_ray_casting.py
"""

import os
import glob
import re
import shutil
from typing import List, Tuple, Dict, Any


def parse_cells_and_spots(lines: List[str]) -> List[Dict[str, Any]]:
    """
    Parse FISH-QUANT file to extract cells with their outlines and spots.
    Handle cases where some cells may not have SPOTS sections.
    
    Args:
        lines: Lines from the FISH-QUANT file
    
    Returns:
        List of cell dictionaries, each containing outline and spots
    """
    cells = []
    current_cell = None
    spot_section_started = False
    header_columns = []
    
    for line in lines:
        line = line.strip()
        
        # Detect new cell
        if line.startswith('CELL'):
            if current_cell is not None:
                cells.append(current_cell)
            cell_name = line.split('\t')[1] if '\t' in line else line.replace('CELL', '').strip()
            current_cell = {
                'name': cell_name,
                'x_coords': [],
                'y_coords': [],
                'spots': []
            }
            spot_section_started = False
            header_columns = []
            continue
        
        if current_cell is None:
            continue
            
        # When we see another CELL while in spot section, stop spot parsing
        if spot_section_started and line.startswith('CELL'):
            spot_section_started = False
            # This line will be processed in the next iteration
            
        # Parse outline coordinates for current cell
        if line.startswith('X_POS') and not spot_section_started:
            coords_str = line.split('\t', 1)[1]
            coords = [float(x) for x in coords_str.split('\t') if x != 'END']
            current_cell['x_coords'] = coords
        elif line.startswith('Y_POS') and not spot_section_started:
            coords_str = line.split('\t', 1)[1]
            coords = [float(y) for y in coords_str.split('\t') if y != 'END']
            current_cell['y_coords'] = coords
        
        # Handle SPOTS section
        elif line == 'SPOTS':
            spot_section_started = True
            header_columns = []
            continue
        
        if spot_section_started:
            # Check for end of spot section (new cell or new coordinate section)
            if line.startswith(('CELL', 'X_POS', 'Y_POS')):
                spot_section_started = False
                continue
                
            if not header_columns:
                # Check if this line is the header or first data line
                try:
                    # Try to parse first column as float to detect data line
                    float(line.split('\t')[0])
                    # If successful, this is data without header - create default header
                    header_columns = [f'col_{i}' for i in range(len(line.split('\t')))]
                    # Process this line as data (don't continue)
                except ValueError:
                    # This is a header line
                    header_columns = line.split('\t')
                    continue
            
            # Parse spot data
            if (line and not line.startswith('#') and 
                not line.startswith('X_POS') and 
                not line.startswith('Y_POS') and
                not line.startswith('CELL')):
                
                values = line.split('\t')
                if len(values) >= len(header_columns):
                    try:
                        float(values[0])  # Test if first column is numeric
                        spot_dict = {}
                        for i, col in enumerate(header_columns):
                            try:
                                # Keep original string format for large numbers to preserve scientific notation
                                val_str = values[i] if i < len(values) else ''
                                val_float = float(val_str)
                                # If it's in scientific notation format or very large, keep as string
                                if 'e' in val_str.lower() or abs(val_float) >= 1e6:
                                    spot_dict[col] = val_str
                                else:
                                    spot_dict[col] = val_float
                            except (ValueError, IndexError):
                                spot_dict[col] = values[i] if i < len(values) else ''
                        current_cell['spots'].append(spot_dict)
                    except ValueError:
                        continue
    
    # Add the last cell
    if current_cell is not None:
        cells.append(current_cell)
    
    return cells


def point_in_polygon_ray_casting(point: Tuple[float, float], polygon: List[Tuple[float, float]]) -> bool:
    """
    Determine if a point is inside a polygon using the ray casting algorithm.
    
    Args:
        point: (x, y) coordinates of the test point
        polygon: List of (x, y) coordinates defining the polygon vertices
    
    Returns:
        True if point is inside polygon, False otherwise
    """
    x, y = point
    n = len(polygon)
    inside = False
    
    p1x, p1y = polygon[0]
    for i in range(1, n + 1):
        p2x, p2y = polygon[i % n]
        
        # Check if point is on the same horizontal level as the edge
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    # Calculate intersection if edge is not horizontal
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    # Check if point is to the left of intersection
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        
        p1x, p1y = p2x, p2y
    
    return inside


def write_filtered_file_multicell(original_lines: List[str], cells: List[Dict[str, Any]], 
                                 all_filtered_spots: List[Dict[str, Any]], output_file: str):
    """
    Write filtered spots to output file, preserving original multi-cell format.
    Copy original file line by line, only replacing spot data sections.
    """
    # Create a mapping from cell name to filtered spots
    cell_spots_map = {}
    for cell in cells:
        cell_spots_map[cell['name']] = cell.get('filtered_spots', [])
    
    with open(output_file, 'w', encoding='utf-8') as f:
        current_cell_name = None
        in_spots_section = False
        header_columns = []
        spot_data_started = False
        
        i = 0
        while i < len(original_lines):
            line = original_lines[i]
            line_stripped = line.strip()
            
            # Detect CELL lines
            if line_stripped.startswith('CELL'):
                current_cell_name = line_stripped.split('\t')[1] if '\t' in line_stripped else line_stripped.replace('CELL', '').strip()
                in_spots_section = False
                spot_data_started = False
                header_columns = []
                f.write(line)
                i += 1
                continue
            
            # Detect SPOTS section
            if line_stripped == 'SPOTS':
                in_spots_section = True
                spot_data_started = False
                header_columns = []
                f.write(line)
                i += 1
                continue
            
            # Handle lines in SPOTS section
            if in_spots_section:
                # Check if this might be the header line
                if not spot_data_started and not header_columns and line_stripped:
                    try:
                        # Test if first column is numeric (spot data)
                        float(line_stripped.split('\t')[0])
                        # It's data, so no explicit header was provided
                        header_columns = [f'col_{j}' for j in range(len(line_stripped.split('\t')))]
                        spot_data_started = True
                        # Process this line as data (don't write original, write filtered)
                    except ValueError:
                        # It's a header line
                        header_columns = line_stripped.split('\t')
                        f.write(line)
                        i += 1
                        continue
                
                # Handle spot data lines
                if spot_data_started or (header_columns and line_stripped):
                    # Check if this is spot data
                    try:
                        float(line_stripped.split('\t')[0])
                        # This is spot data - replace with filtered data
                        if not spot_data_started:
                            spot_data_started = True
                        
                        # Write filtered spots for current cell (only once)
                        if current_cell_name in cell_spots_map:
                            filtered_spots = cell_spots_map[current_cell_name]
                            for spot in filtered_spots:
                                row_values = []
                                for col in header_columns:
                                    value = spot.get(col, '')
                                    if isinstance(value, float):
                                        if value == int(value):
                                            row_values.append(str(int(value)))
                                        elif abs(value) >= 1e6 or (abs(value) < 1e-3 and value != 0):
                                            row_values.append(f"{value:.5e}")
                                        else:
                                            row_values.append(f"{value:.6g}")
                                    elif isinstance(value, str) and value:
                                        row_values.append(value)
                                    else:
                                        row_values.append(str(value))
                                f.write('\t'.join(row_values) + '\n')
                            # Remove from map so we don't write again
                            del cell_spots_map[current_cell_name]
                        
                        # Skip original spot data until next section
                        while i < len(original_lines):
                            line = original_lines[i]
                            line_stripped = line.strip()
                            if line_stripped.startswith(('CELL', 'X_POS', 'Y_POS', 'SPOTS')) or not line_stripped:
                                in_spots_section = False
                                break
                            # Skip spot data lines
                            try:
                                float(line_stripped.split('\t')[0])
                                i += 1
                                continue
                            except (ValueError, IndexError):
                                break
                        continue
                        
                    except (ValueError, IndexError):
                        # Not spot data, write as-is
                        f.write(line)
                        i += 1
                        continue
            
            # For all other lines, write as-is
            f.write(line)
            i += 1


def get_pixel_size(lines: List[str]) -> Tuple[float, float]:
    """
    Extract pixel size from FISH-QUANT file.
    
    Args:
        lines: Lines from the FISH-QUANT file
    
    Returns:
        Tuple of (pixel_xy_nm, pixel_z_nm)
    """
    for i, line in enumerate(lines):
        if line.strip() == 'PARAMETERS':
            # Next line should contain the headers
            if i + 1 < len(lines):
                headers = lines[i + 1].strip().split('\t')
                if i + 2 < len(lines):
                    values = lines[i + 2].strip().split('\t')
                    
                    # Find Pix-XY and Pix-Z columns
                    pixel_xy_nm = 64.5  # Default value
                    pixel_z_nm = 200    # Default value
                    
                    for j, header in enumerate(headers):
                        if header == 'Pix-XY' and j < len(values):
                            pixel_xy_nm = float(values[j])
                        elif header == 'Pix-Z' and j < len(values):
                            pixel_z_nm = float(values[j])
                    
                    return pixel_xy_nm, pixel_z_nm
    
    return 64.5, 200  # Default values


def filter_spots_in_file(input_file: str, output_file: str) -> Dict[str, int]:
    """
    Filter spots in a single FISH-QUANT file.
    
    Args:
        input_file: Path to input FISH-QUANT file
        output_file: Path to output filtered file
    
    Returns:
        Dictionary with statistics (total_spots, filtered_spots, removed_spots)
    """
    # Read the file
    with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    # Keep original lines with their newline characters
    # lines = [line.rstrip('\n\r') for line in lines]
    
    # Extract pixel size
    pixel_xy_nm, pixel_z_nm = get_pixel_size(lines)
    
    # Parse cells with their outlines and spots
    cells = parse_cells_and_spots(lines)
    
    if not cells:
        print(f"Warning: No cells found in {input_file}")
        shutil.copy2(input_file, output_file)
        return {'total_spots': 0, 'filtered_spots': 0, 'removed_spots': 0}
    
    total_spots = 0
    all_filtered_spots = []
    
    # Process each cell
    for cell in cells:
        x_coords = cell['x_coords']
        y_coords = cell['y_coords']
        spots = cell['spots']
        
        if not x_coords or not y_coords:
            print(f"  Warning: No outline coordinates for cell {cell['name']}")
            # Include all spots if no outline
            cell['filtered_spots'] = spots
            all_filtered_spots.extend(spots)
            total_spots += len(spots)
            continue
        
        # Create polygon from outline coordinates
        polygon_px = list(zip(x_coords, y_coords))
        
        cell_total = len(spots)
        cell_filtered = []
        
        # Filter spots using ray casting
        for spot in spots:
            try:
                # Convert spot coordinates from nm to pixels
                pos_x = spot.get('Pos_X', 0)
                pos_y = spot.get('Pos_Y', 0)
                
                # Ensure coordinates are numeric
                if isinstance(pos_x, str):
                    pos_x = float(pos_x)
                if isinstance(pos_y, str):
                    pos_y = float(pos_y)
                
                spot_x_px = pos_x / pixel_xy_nm
                spot_y_px = pos_y / pixel_xy_nm
                
                # Check if spot is inside polygon
                if point_in_polygon_ray_casting((spot_x_px, spot_y_px), polygon_px):
                    cell_filtered.append(spot)
                    
            except (ValueError, TypeError, ZeroDivisionError) as e:
                print(f"  Warning: Skipping invalid spot data: {e}")
                continue
        
        # Store filtered spots for this cell
        cell['filtered_spots'] = cell_filtered
        all_filtered_spots.extend(cell_filtered)
        total_spots += cell_total
        
        print(f"  Cell {cell['name']}: {len(cell_filtered)}/{cell_total} spots kept")
    
    filtered_count = len(all_filtered_spots)
    removed_count = total_spots - filtered_count
    
    # Write filtered results to output file
    write_filtered_file_multicell(lines, cells, all_filtered_spots, output_file)
    
    return {
        'total_spots': total_spots,
        'filtered_spots': filtered_count,
        'removed_spots': removed_count
    }





def main():
    """Main function to process all spot files."""
    # Define input and output directories
    input_dir = "/Volumes/ExFAT/mRNA_molecules/Y333 ATP6 ATP2/atp6_corrected_spots"
    output_dir = "/Volumes/ExFAT/mRNA_molecules/Y333 ATP6 ATP2/atp6_filtered_spots"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all .txt spot files (excluding system files and settings files)
    spot_files = glob.glob(os.path.join(input_dir, "*.txt"))
    spot_files = [f for f in spot_files if not os.path.basename(f).startswith('._') 
                  and 'batch_settings' not in os.path.basename(f)]
    
    if not spot_files:
        print(f"No .txt files found in {input_dir}")
        return
    
    print(f"Found {len(spot_files)} spot files to process")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    total_stats = {'total_spots': 0, 'filtered_spots': 0, 'removed_spots': 0}
    
    # Process each file
    for i, input_file in enumerate(sorted(spot_files), 1):
        filename = os.path.basename(input_file)
        output_file = os.path.join(output_dir, filename)
        
        print(f"Processing [{i}/{len(spot_files)}]: {filename}")
        
        try:
            stats = filter_spots_in_file(input_file, output_file)
            
            print(f"  Total spots: {stats['total_spots']}")
            print(f"  Filtered spots: {stats['filtered_spots']}")
            print(f"  Removed spots: {stats['removed_spots']}")
            if stats['total_spots'] > 0:
                removal_percent = (stats['removed_spots'] / stats['total_spots']) * 100
                print(f"  Removal rate: {removal_percent:.1f}%")
            print()
            
            # Update total statistics
            for key in total_stats:
                total_stats[key] += stats[key]
            
        except Exception as e:
            print(f"  Error processing {filename}: {str(e)}")
            print()
    
    # Print summary statistics
    print("="*50)
    print("SUMMARY STATISTICS")
    print("="*50)
    print(f"Files processed: {len(spot_files)}")
    print(f"Total spots: {total_stats['total_spots']}")
    print(f"Filtered spots: {total_stats['filtered_spots']}")
    print(f"Removed spots: {total_stats['removed_spots']}")
    if total_stats['total_spots'] > 0:
        overall_removal_percent = (total_stats['removed_spots'] / total_stats['total_spots']) * 100
        print(f"Overall removal rate: {overall_removal_percent:.1f}%")
    print(f"\nFiltered files saved to: {output_dir}")


if __name__ == "__main__":
    main() 