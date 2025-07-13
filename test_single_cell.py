#!/usr/bin/env python3
"""
Test single cell analysis to debug issues
"""

import sys
import os
sys.path.append('.')

from batch_analyze_alignment import find_all_cell_files, analyze_single_cell

def test_single_cell():
    """Test analysis of a single cell"""
    print("Finding all cell files...")
    cell_files = find_all_cell_files()
    
    if not cell_files:
        print("No cell files found!")
        return
    
    print(f"Found {len(cell_files)} cells")
    print("Testing first cell...")
    
    # Test the first cell
    first_cell = cell_files[0]
    print(f"Testing: {first_cell['experiment']}_{first_cell['position']}_cell_{first_cell['cell_id']}")
    print(f"Skeleton file: {first_cell['skeleton_file']}")
    print(f"Spots file: {first_cell['spots_file']}")
    print(f"Mapping file: {first_cell['mapping_file']}")
    
    # Test with quiet=False to see all output
    try:
        result = analyze_single_cell(first_cell, method='static', quiet=False)
        if result:
            print(f"SUCCESS! Result keys: {list(result.keys())}")
            print(f"Number of spots: {result['n_spots']}")
            print(f"Number of skeleton points: {result['n_skeleton_points']}")
            print(f"Mean distance: {result['mean_distance']:.3f}")
        else:
            print("FAILED: analyze_single_cell returned None")
    except Exception as e:
        print(f"ERROR in analyze_single_cell: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_single_cell() 