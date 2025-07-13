#!/usr/bin/env python3
"""
Script to replace slices in a multi-slice image with blank ones.
Keeps only a specified range of slices and sets all others to blank/zero.

Usage:
    # Single file processing
    python replace_slices.py input_image.tif output_image.tif --keep-range 10 20
    python replace_slices.py input_image.tif output_image.tif --keep-slices 5,10,15,20
    
    # Batch processing
    python replace_slices.py --batch input_folder output_folder --keep-range 10 20
    python replace_slices.py --batch input_folder output_folder --keep-slices 5,10,15,20 --pattern "*.tif"
"""

import argparse
import numpy as np
from skimage import io
from pathlib import Path
import sys
import glob
import os
from tqdm import tqdm

def load_multi_slice_image(filepath):
    """Load a multi-slice image from file."""
    try:
        # Try loading with skimage first (handles TIFF stacks well)
        image = io.imread(filepath)
        print(f"Loaded image with shape: {image.shape}")
        return image
    except Exception as e:
        print(f"Error loading image {filepath}: {e}")
        return None

def save_multi_slice_image(image, filepath):
    """Save a multi-slice image to file."""
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        io.imsave(filepath, image)
        print(f"Saved image to: {filepath}")
        return True
    except Exception as e:
        print(f"Error saving image {filepath}: {e}")
        return False

def replace_slices_with_blank(image, keep_slices):
    """
    Replace all slices except the specified ones with blank (zero) slices.
    
    Args:
        image: numpy array of shape (slices, height, width) or (slices, height, width, channels)
        keep_slices: list of slice indices to keep (0-based)
    
    Returns:
        Modified image with non-selected slices set to zero
    """
    if len(image.shape) < 3:
        raise ValueError("Image must have at least 3 dimensions (slices, height, width)")
    
    total_slices = image.shape[0]
    
    # Validate slice indices
    valid_slices = [s for s in keep_slices if 0 <= s < total_slices]
    invalid_slices = [s for s in keep_slices if s not in valid_slices]
    
    if invalid_slices:
        print(f"Warning: Invalid slice indices (out of range): {invalid_slices}")
    
    if not valid_slices:
        raise ValueError("No valid slice indices provided")
    
    # Create a copy of the image
    result = image.copy()
    
    # Set all slices to blank first
    result[:] = 0
    
    # Restore the slices we want to keep
    for slice_idx in valid_slices:
        result[slice_idx] = image[slice_idx]
    
    return result

def find_image_files(input_folder, pattern="*"):
    """Find all image files in the input folder matching the pattern."""
    input_path = Path(input_folder)
    if not input_path.exists():
        raise ValueError(f"Input folder does not exist: {input_folder}")
    
    # Common image extensions (both upper and lowercase)
    image_extensions = ['.tif', '.tiff', '.png', '.jpg', '.jpeg', '.bmp', '.npy',
                       '.TIF', '.TIFF', '.PNG', '.JPG', '.JPEG', '.BMP', '.NPY']
    
    image_files = []
    
    # If pattern is specific (contains an extension), use it directly
    if any(ext in pattern for ext in image_extensions):
        files = list(input_path.glob(pattern))
        image_files.extend(files)
    else:
        # If pattern is generic (like "*"), try with all extensions
        for ext in image_extensions:
            search_pattern = pattern.replace('*', f'*{ext}') if '*' in pattern else f'{pattern}{ext}'
            files = list(input_path.glob(search_pattern))
            image_files.extend(files)
    
    # Filter out macOS metadata files and other unwanted files
    image_files = [f for f in image_files if not f.name.startswith('._')]
    
    # Remove duplicates and sort
    image_files = sorted(list(set(image_files)))
    
    return image_files

def process_single_image(input_file, output_file, keep_slices, verbose=True):
    """Process a single image file."""
    if verbose:
        print(f"\nProcessing: {input_file}")
    
    # Load the image
    image = load_multi_slice_image(input_file)
    if image is None:
        return False
    
    try:
        # Process the image
        result = replace_slices_with_blank(image, keep_slices)
        
        # Save the result
        success = save_multi_slice_image(result, output_file)
        
        if verbose and success:
            total_slices = image.shape[0]
            kept_slices = len([s for s in keep_slices if 0 <= s < total_slices])
            print(f"Kept {kept_slices}/{total_slices} slices")
        
        return success
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return False

def batch_process(input_folder, output_folder, keep_slices, pattern="*", suffix="_processed"):
    """Process all images in a folder."""
    print(f"Batch processing folder: {input_folder}")
    print(f"Output folder: {output_folder}")
    print(f"File pattern: {pattern}")
    
    # Find all image files
    image_files = find_image_files(input_folder, pattern)
    
    if not image_files:
        print(f"No image files found in {input_folder} matching pattern '{pattern}'")
        return
    
    print(f"Found {len(image_files)} image files to process")
    
    # Create output folder
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Process each file with progress bar
    successful = 0
    failed = 0
    
    for input_file in tqdm(image_files, desc="Processing images"):
        # Generate output filename
        input_path = Path(input_file)
        output_filename = input_path.stem + suffix + input_path.suffix
        output_file = output_path / output_filename
        
        # Process the image
        success = process_single_image(input_file, output_file, keep_slices, verbose=False)
        
        if success:
            successful += 1
        else:
            failed += 1
            print(f"Failed to process: {input_file}")
    
    print(f"\nBatch processing complete:")
    print(f"Successfully processed: {successful}")
    print(f"Failed: {failed}")

def parse_slice_list(slice_string):
    """Parse comma-separated slice indices."""
    try:
        return [int(s.strip()) for s in slice_string.split(',')]
    except ValueError as e:
        raise ValueError(f"Invalid slice list format: {slice_string}. Use comma-separated integers.")

def main():
    parser = argparse.ArgumentParser(
        description="Replace slices in multi-slice image with blank ones (single file or batch processing)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file processing
  python replace_slices.py input.tif output.tif --keep-range 10 20
  python replace_slices.py input.tif output.tif --keep-slices 5,10,15,20
  
  # Batch processing - process all images in folder
  python replace_slices.py --batch input_folder output_folder --keep-range 10 20
  python replace_slices.py --batch input_folder output_folder --keep-slices 5,10,15,20
  
  # Batch processing with file pattern and custom suffix
  python replace_slices.py --batch input_folder output_folder --keep-range 10 20 --pattern "*.tif" --suffix "_blank"
  
  # Preview mode (works for both single and batch)
  python replace_slices.py input.tif output.tif --keep-range 10 20 --preview
  python replace_slices.py --batch input_folder output_folder --keep-range 10 20 --preview
        """
    )
    
    # Batch processing option
    parser.add_argument("--batch", action="store_true",
                       help="Enable batch processing mode")
    
    # Input/output arguments - different behavior for batch vs single
    parser.add_argument("input", help="Input file (single mode) or input folder (batch mode)")
    parser.add_argument("output", help="Output file (single mode) or output folder (batch mode)")
    
    # Slice selection (mutually exclusive)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--keep-range", nargs=2, type=int, metavar=("START", "END"),
                      help="Keep slice range (inclusive). E.g., --keep-range 10 20")
    group.add_argument("--keep-slices", type=str, metavar="SLICE_LIST",
                      help="Keep specific slices (comma-separated). E.g., --keep-slices 5,10,15,20")
    
    # Batch processing options
    parser.add_argument("--pattern", default="*", 
                       help="File pattern for batch processing (default: '*' for all image files)")
    parser.add_argument("--suffix", default="_processed",
                       help="Suffix to add to output filenames in batch mode (default: '_processed')")
    
    # General options
    parser.add_argument("--preview", action="store_true",
                       help="Preview which slices will be kept without saving")
    
    args = parser.parse_args()
    
    # Determine which slices to keep
    if args.keep_range:
        start, end = args.keep_range
        if start > end:
            print("Error: Start slice must be <= end slice")
            sys.exit(1)
        keep_slices = list(range(start, end + 1))  # inclusive range
        print(f"Keeping slice range: {start} to {end} (inclusive)")
    else:
        keep_slices = parse_slice_list(args.keep_slices)
        print(f"Keeping specific slices: {sorted(keep_slices)}")
    
    if args.batch:
        # Batch processing mode
        if not Path(args.input).exists():
            print(f"Error: Input folder {args.input} does not exist")
            sys.exit(1)
        
        if args.preview:
            # Preview batch processing
            image_files = find_image_files(args.input, args.pattern)
            print(f"\nBatch Preview Mode:")
            print(f"Input folder: {args.input}")
            print(f"Output folder: {args.output}")
            print(f"File pattern: {args.pattern}")
            print(f"Files to process: {len(image_files)}")
            print(f"Slices to keep: {sorted(keep_slices)}")
            print(f"Output suffix: {args.suffix}")
            if image_files:
                print("Files found:")
                for f in image_files[:10]:  # Show first 10 files
                    print(f"  {f}")
                if len(image_files) > 10:
                    print(f"  ... and {len(image_files) - 10} more files")
        else:
            batch_process(args.input, args.output, keep_slices, args.pattern, args.suffix)
    
    else:
        # Single file processing mode
        if not Path(args.input).exists():
            print(f"Error: Input file {args.input} does not exist")
            sys.exit(1)
        
        # Load the image for preview or processing
        print(f"Loading image: {args.input}")
        image = load_multi_slice_image(args.input)
        if image is None:
            sys.exit(1)
        
        if args.preview:
            total_slices = image.shape[0]
            print(f"\nPreview mode:")
            print(f"Input file: {args.input}")
            print(f"Output file: {args.output}")
            print(f"Total slices: {total_slices}")
            print(f"Slices to keep: {sorted(keep_slices)}")
            print(f"Slices to blank: {total_slices - len(keep_slices)}")
            blank_slices = [i for i in range(total_slices) if i not in keep_slices]
            if len(blank_slices) <= 20:  # Don't print too many
                print(f"Blank slice indices: {blank_slices}")
            else:
                print(f"Blank slice indices: {blank_slices[:10]}...{blank_slices[-10:]} (showing first/last 10)")
        else:
            # Process single file
            success = process_single_image(args.input, args.output, keep_slices)
            if success:
                print("Done!")
            else:
                sys.exit(1)

if __name__ == "__main__":
    main() 