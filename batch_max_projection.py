#!/usr/bin/env python3
"""
Batch Maximum Projection Script with Brightness/Contrast Reset

This script performs maximum projection on image stacks from red and green channel folders
and creates merged images with automatic brightness and contrast normalization. 
Both channels are automatically converted to 16-bit format before merging to ensure 
consistent bit depth. It excludes files with 'filtered_batch' postfix.

Features:
    - Converts all images to 16-bit format before processing
    - Maximum intensity projection along Z-axis
    - Brightness/contrast normalization with multiple methods
    - Creates merged RGB images (red + green channels)
    - Saves individual normalized projections

Usage:
    python batch_max_projection.py <red_channel_folder> <green_channel_folder> [output_folder] [--normalize METHOD]

Normalization Methods:
    - auto_contrast: Automatic contrast adjustment using percentiles (default)
    - full_range: Stretch to full 0-65535 range  
    - histogram_eq: Global histogram equalization
    - adaptive_eq: Adaptive histogram equalization (CLAHE)

Examples:
    python batch_max_projection.py /path/to/red /path/to/green
    python batch_max_projection.py /path/to/red /path/to/green /path/to/output --normalize auto_contrast
    python batch_max_projection.py /path/to/red /path/to/green --normalize full_range
"""

import os
import sys
import glob
import numpy as np
from pathlib import Path
import argparse
from tifffile import imread, imwrite
from skimage import exposure
import re

def find_matching_files(red_folder, green_folder, exclude_pattern="filtered_batch"):
    """
    Find matching red and green channel files based on filename patterns.
    Excludes files containing the exclude_pattern.
    
    Expected patterns:
    Red: yWL333_cy3_ATP2_cy5_ATP6MS2_X_w1CY5-100-_sY.TIF
    Green: deconv_yWL333_cy3_ATP2_cy5_ATP6MS2_X_w4FITC-100-_sY.TIF
    """
    # Get all TIF files from both folders
    red_files = glob.glob(os.path.join(red_folder, "*.TIF")) + glob.glob(os.path.join(red_folder, "*.tif"))
    green_files = glob.glob(os.path.join(green_folder, "*.TIF")) + glob.glob(os.path.join(green_folder, "*.tif"))
    
    # Filter out files with excluded pattern
    red_files = [f for f in red_files if exclude_pattern not in os.path.basename(f)]
    green_files = [f for f in green_files if exclude_pattern not in os.path.basename(f)]
    
    print(f"Found {len(red_files)} red channel files")
    print(f"Found {len(green_files)} green channel files")
    
    # Extract base names for matching
    red_bases = {}
    green_bases = {}
    
    for file in red_files:
        basename = os.path.basename(file)
        # Extract core identifier: experiment number and slide number
        # Pattern: yWL333_cy3_ATP2_cy5_ATP6MS2_X_w1CY5-100-_sY.TIF -> yWL333_cy3_ATP2_cy5_ATP6MS2_X_sY
        match = re.search(r'(yWL333_cy3_ATP2_cy5_ATP6MS2_\d+)_w\d+[A-Z0-9-]+_(s\d+)', basename)
        if match:
            core_name = f"{match.group(1)}_{match.group(2)}"
            red_bases[core_name] = file
            print(f"Red: {basename} -> {core_name}")
    
    for file in green_files:
        basename = os.path.basename(file)
        # Extract core identifier from green files (remove deconv_ prefix)
        # Pattern: deconv_yWL333_cy3_ATP2_cy5_ATP6MS2_X_w4FITC-100-_sY.TIF -> yWL333_cy3_ATP2_cy5_ATP6MS2_X_sY
        basename_no_prefix = re.sub(r'^deconv_', '', basename)
        match = re.search(r'(yWL333_cy3_ATP2_cy5_ATP6MS2_\d+)_w\d+[A-Z0-9-]+_(s\d+)', basename_no_prefix)
        if match:
            core_name = f"{match.group(1)}_{match.group(2)}"
            green_bases[core_name] = file
            print(f"Green: {basename} -> {core_name}")
    
    # Find matching pairs
    matches = []
    for core_name in red_bases:
        if core_name in green_bases:
            matches.append((red_bases[core_name], green_bases[core_name]))
            print(f"Match found: {core_name}")
        else:
            print(f"No green match for red: {core_name}")
    
    # Check for unmatched green files
    for core_name in green_bases:
        if core_name not in red_bases:
            print(f"No red match for green: {core_name}")
    
    return matches

def max_projection(image_stack):
    """
    Perform maximum intensity projection along the Z-axis.
    """
    if image_stack.ndim == 3:
        return np.max(image_stack, axis=0)
    elif image_stack.ndim == 2:
        return image_stack
    else:
        raise ValueError(f"Unexpected image dimensions: {image_stack.shape}")

def convert_to_16bit(image):
    """
    Convert image to 16-bit format regardless of input bit depth.
    """
    # Determine the current bit depth
    if image.dtype == np.uint8:
        # Convert 8-bit to 16-bit by scaling
        return (image.astype(np.float32) * 257).astype(np.uint16)  # 257 = 65535/255
    elif image.dtype == np.uint16:
        # Already 16-bit
        return image
    elif image.dtype == np.uint32:
        # Convert 32-bit to 16-bit by scaling down
        return (image.astype(np.float64) / 65537).astype(np.uint16)  # 65537 = 2^32/2^16
    elif image.dtype in [np.float32, np.float64]:
        # Convert float to 16-bit
        # Assume float range is 0-1, scale to 0-65535
        normalized = np.clip(image, 0, 1)
        return (normalized * 65535).astype(np.uint16)
    else:
        # Force conversion to 16-bit
        return image.astype(np.uint16)

def normalize_image(image, method='auto_contrast', percentile_range=(1, 99)):
    """
    Normalize image intensity with various brightness and contrast adjustment methods.
    Ensures output is 16-bit format.
    
    Methods:
    - 'auto_contrast': Automatic contrast adjustment using percentiles
    - 'full_range': Stretch to full 0-65535 range
    - 'histogram_eq': Histogram equalization
    - 'adaptive_eq': Adaptive histogram equalization
    """
    
    # First convert to 16-bit if not already
    image_16bit = convert_to_16bit(image)
    
    if method == 'auto_contrast':
        # Auto contrast using percentiles to reset brightness/contrast
        p_low, p_high = np.percentile(image_16bit, percentile_range)
        if p_high > p_low:  # Avoid division by zero
            normalized = exposure.rescale_intensity(
                image_16bit, 
                in_range=(p_low, p_high), 
                out_range=(0, 65535)
            ).astype(np.uint16)
        else:
            normalized = image_16bit
            
    elif method == 'full_range':
        # Stretch to full dynamic range
        min_val, max_val = image_16bit.min(), image_16bit.max()
        if max_val > min_val:
            normalized = exposure.rescale_intensity(
                image_16bit, 
                in_range=(min_val, max_val), 
                out_range=(0, 65535)
            ).astype(np.uint16)
        else:
            normalized = image_16bit
            
    elif method == 'histogram_eq':
        # Global histogram equalization
        # Convert to float for processing, then back to 16-bit
        image_float = image_16bit.astype(np.float64) / 65535.0
        normalized = exposure.equalize_hist(image_float)
        normalized = (normalized * 65535).astype(np.uint16)
        
    elif method == 'adaptive_eq':
        # Adaptive histogram equalization (CLAHE)
        # Convert to float for processing, then back to 16-bit
        image_float = image_16bit.astype(np.float64) / 65535.0
        normalized = exposure.equalize_adapthist(image_float, clip_limit=0.02)
        normalized = (normalized * 65535).astype(np.uint16)
        
    else:
        # Default to auto contrast
        p_low, p_high = np.percentile(image_16bit, percentile_range)
        if p_high > p_low:
            normalized = exposure.rescale_intensity(
                image_16bit, 
                in_range=(p_low, p_high), 
                out_range=(0, 65535)
            ).astype(np.uint16)
        else:
            normalized = image_16bit
    
    return normalized

def create_merged_image(red_proj, green_proj, normalization_method='auto_contrast'):
    """
    Create a merged RGB image from red and green projections with brightness/contrast reset.
    Both input channels must be 16-bit format. Output is 16-bit RGB.
    """
    # Verify both inputs are 16-bit
    assert red_proj.dtype == np.uint16, f"Red channel must be 16-bit, got {red_proj.dtype}"
    assert green_proj.dtype == np.uint16, f"Green channel must be 16-bit, got {green_proj.dtype}"
    
    # Normalize both channels with specified method (input already 16-bit)
    red_norm = normalize_image(red_proj, method=normalization_method)
    green_norm = normalize_image(green_proj, method=normalization_method)
    
    # Verify normalization preserved 16-bit format
    assert red_norm.dtype == np.uint16, f"Red normalization failed to maintain 16-bit"
    assert green_norm.dtype == np.uint16, f"Green normalization failed to maintain 16-bit"
    
    # Create 16-bit RGB image
    height, width = red_norm.shape
    merged = np.zeros((height, width, 3), dtype=np.uint16)
    
    # Assign channels: Red=0, Green=1, Blue=2
    merged[:, :, 0] = red_norm    # Red channel
    merged[:, :, 1] = green_norm  # Green channel
    merged[:, :, 2] = 0           # Blue channel (empty)
    
    return merged

def process_image_pair(red_file, green_file, output_folder, base_name, normalization_method='auto_contrast'):
    """
    Process a pair of red and green channel images with brightness/contrast reset.
    Ensures both channels are converted to 16-bit before merging.
    """
    try:
        print(f"Processing: {base_name}")
        
        # Read images
        red_stack = imread(red_file)
        green_stack = imread(green_file)
        
        print(f"  Red shape: {red_stack.shape}, dtype: {red_stack.dtype}")
        print(f"  Green shape: {green_stack.shape}, dtype: {green_stack.dtype}")
        
        # Perform maximum projection
        red_proj = max_projection(red_stack)
        green_proj = max_projection(green_stack)
        
        # Convert both projections to 16-bit before processing
        red_proj_16bit = convert_to_16bit(red_proj)
        green_proj_16bit = convert_to_16bit(green_proj)
        
        print(f"  Converted to 16-bit - Red: {red_proj_16bit.dtype}, Green: {green_proj_16bit.dtype}")
        
        # Create merged image with specified normalization (both channels now 16-bit)
        merged = create_merged_image(red_proj_16bit, green_proj_16bit, normalization_method)
        
        # Save individual projections with brightness/contrast reset
        red_output = os.path.join(output_folder, f"{base_name}_red_maxproj.tif")
        green_output = os.path.join(output_folder, f"{base_name}_green_maxproj.tif")
        merged_output = os.path.join(output_folder, f"{base_name}_merged_maxproj.tif")
        
        # Normalize and save (input is already 16-bit)
        red_normalized = normalize_image(red_proj_16bit, method=normalization_method)
        green_normalized = normalize_image(green_proj_16bit, method=normalization_method)
        
        imwrite(red_output, red_normalized)
        imwrite(green_output, green_normalized)
        imwrite(merged_output, merged)
        
        print(f"  Saved: {os.path.basename(merged_output)} (16-bit)")
        print(f"  Normalization: {normalization_method}")
        print(f"  Output bit depth - Red: {red_normalized.dtype}, Green: {green_normalized.dtype}, Merged: {merged.dtype}")
        
    except Exception as e:
        print(f"Error processing {base_name}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Batch Maximum Projection for Red and Green Channels with Brightness/Contrast Reset')
    parser.add_argument('red_folder', help='Path to red channel folder')
    parser.add_argument('green_folder', help='Path to green channel folder')
    parser.add_argument('output_folder', nargs='?', default='max_projections', 
                       help='Output folder (default: max_projections)')
    parser.add_argument('--exclude', default='filtered_batch', 
                       help='Pattern to exclude from processing (default: filtered_batch)')
    parser.add_argument('--normalize', default='auto_contrast', 
                       choices=['auto_contrast', 'full_range', 'histogram_eq', 'adaptive_eq'],
                       help='Brightness/contrast normalization method (default: auto_contrast)')
    
    args = parser.parse_args()
    
    # Validate input folders
    if not os.path.isdir(args.red_folder):
        print(f"Error: Red channel folder '{args.red_folder}' does not exist.")
        sys.exit(1)
    
    if not os.path.isdir(args.green_folder):
        print(f"Error: Green channel folder '{args.green_folder}' does not exist.")
        sys.exit(1)
    
    # Create output folder
    os.makedirs(args.output_folder, exist_ok=True)
    print(f"Output folder: {args.output_folder}")
    
    # Find matching files
    print("Finding matching red and green channel files...")
    matches = find_matching_files(args.red_folder, args.green_folder, args.exclude)
    
    if not matches:
        print("No matching red/green channel pairs found!")
        print(f"Red folder: {args.red_folder}")
        print(f"Green folder: {args.green_folder}")
        print(f"Excluding pattern: {args.exclude}")
        sys.exit(1)
    
    print(f"Found {len(matches)} matching pairs")
    print(f"Using normalization method: {args.normalize}")
    
    # Process each matching pair
    for i, (red_file, green_file) in enumerate(matches, 1):
        # Generate base name for output
        red_base = os.path.splitext(os.path.basename(red_file))[0]
        base_name = re.sub(r'_w\d+[A-Z0-9-]+', '', red_base)
        base_name = re.sub(r'^deconv_', '', base_name)
        
        print(f"\n[{i}/{len(matches)}]")
        process_image_pair(red_file, green_file, args.output_folder, base_name, args.normalize)
    
    print(f"\nProcessing complete! Results saved to: {args.output_folder}")
    print(f"All images processed with {args.normalize} brightness/contrast normalization")

if __name__ == "__main__":
    main() 