#!/usr/bin/env python3
"""
Script to fit ellipses to individual cell mask instances from a TIFF image.

This script loads a TIFF image containing cell masks, identifies individual cell instances,
and fits an ellipse to each cell mask. The results are saved to a CSV file with ellipse parameters.
"""

import numpy as np
import pandas as pd
import tifffile
import cv2
from skimage import measure, morphology
from skimage.segmentation import clear_border
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import glob
import os


def load_mask_image(image_path):
    """Load the TIFF mask image."""
    try:
        image = tifffile.imread(image_path)
        print(f"Loaded image with shape: {image.shape}, dtype: {image.dtype}")
        return image
    except Exception as e:
        print(f"Error loading image: {e}")
        return None


def preprocess_mask(mask):
    """
    Preprocess the mask where each cell instance already has a unique value.
    
    Args:
        mask: Input mask image where each cell has a unique pixel value
        
    Returns:
        labeled_mask: Labeled image with each cell having a unique integer label
    """
    # The mask already contains unique values for each cell instance
    # Background should be 0, each cell has a unique non-zero value
    
    # Get unique values (excluding background = 0)
    unique_values = np.unique(mask)
    unique_values = unique_values[unique_values > 0]  # Remove background
    
    print(f"Found {len(unique_values)} cell instances with values: {unique_values[:10]}{'...' if len(unique_values) > 10 else ''}")
    
    # The mask is already properly labeled, just return it
    # But we may want to relabel to ensure consecutive numbering starting from 1
    labeled_mask = np.zeros_like(mask)
    for i, value in enumerate(unique_values, 1):
        labeled_mask[mask == value] = i
    
    print(f"Relabeled to consecutive values 1-{len(unique_values)}")
    
    return labeled_mask


def fit_ellipse_to_mask(mask_region, pixel_size_um=1.0):
    """
    Fit an ellipse to a single cell mask region using OpenCV.
    
    Args:
        mask_region: Binary mask of a single cell
        pixel_size_um: Pixel size in micrometers per pixel
        
    Returns:
        ellipse_params: Dictionary containing ellipse parameters
    """
    # Find contours
    contours, _ = cv2.findContours(
        mask_region.astype(np.uint8), 
        cv2.RETR_EXTERNAL, 
        cv2.CHAIN_APPROX_SIMPLE
    )
    
    if len(contours) == 0:
        return None
    
    # Use the largest contour
    largest_contour = max(contours, key=cv2.contourArea)
    
    # Need at least 5 points to fit an ellipse
    if len(largest_contour) < 5:
        return None
    
    try:
        # Fit ellipse
        ellipse = cv2.fitEllipse(largest_contour)
        
        # Extract parameters
        center_x, center_y = ellipse[0]
        width, height = ellipse[1]  # width and height from OpenCV
        angle = ellipse[2]  # Rotation angle in degrees
        
        # OpenCV returns (width, height) but we need to determine which is major/minor
        # and adjust angle accordingly
        if width > height:
            major_axis = width
            minor_axis = height
            # Angle is already correct for this case
            corrected_angle = angle
        else:
            major_axis = height
            minor_axis = width
            # Need to adjust angle by 90 degrees
            corrected_angle = angle + 90
            if corrected_angle > 180:
                corrected_angle -= 180
        
        # Convert axes to micrometers
        major_axis_um = major_axis * pixel_size_um
        minor_axis_um = minor_axis * pixel_size_um
        
        # Calculate area in µm² and eccentricity
        area_um2 = np.pi * (major_axis_um / 2) * (minor_axis_um / 2)
        eccentricity = np.sqrt(1 - (minor_axis / major_axis)**2)
        
        # Calculate volume in µm³ (ellipsoid rotating around major axis, height = minor axis)
        # V = (4/3) * π * a * b * c, where a = major_axis/2, b = c = minor_axis/2
        volume_um3 = (4/3) * np.pi * (major_axis_um / 2) * (minor_axis_um / 2)**2
        
        return {
            'center_x': center_x * pixel_size_um,
            'center_y': center_y * pixel_size_um,
            'major_axis': major_axis_um,
            'minor_axis': minor_axis_um,
            'angle': corrected_angle,
            'area': area_um2,
            'eccentricity': eccentricity,
            'volume': volume_um3,
            'contour_area': cv2.contourArea(largest_contour) * pixel_size_um**2
        }
        
    except Exception as e:
        print(f"Error fitting ellipse: {e}")
        return None


def fit_ellipses_to_all_cells(labeled_mask, pixel_size_um=1.0):
    """
    Fit ellipses to all cell instances in the labeled mask.
    
    Args:
        labeled_mask: Labeled image with each cell having a unique integer label
        pixel_size_um: Pixel size in micrometers per pixel
        
    Returns:
        results: List of dictionaries containing ellipse parameters for each cell
    """
    results = []
    
    for cell_id in range(1, labeled_mask.max() + 1):
        # Extract individual cell mask
        cell_mask = (labeled_mask == cell_id).astype(np.uint8)
        
        # Fit ellipse
        ellipse_params = fit_ellipse_to_mask(cell_mask, pixel_size_um)
        
        if ellipse_params is not None:
            ellipse_params['cell_id'] = cell_id
            results.append(ellipse_params)
            print(f"Successfully fitted ellipse to cell {cell_id}")
        else:
            print(f"Failed to fit ellipse to cell {cell_id}")
    
    return results


def visualize_ellipses(image, labeled_mask, ellipse_results, output_path, pixel_size_um=1.0):
    """
    Create a visualization showing the original masks and fitted ellipses.
    
    Args:
        image: Original image
        labeled_mask: Labeled mask image
        ellipse_results: List of ellipse parameters
        output_path: Path to save the visualization
        pixel_size_um: Pixel size in micrometers per pixel
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    
    # Show original masks
    axes[0].imshow(labeled_mask, cmap='tab20')
    axes[0].set_title('Original Cell Masks')
    axes[0].axis('off')
    
    # Show ellipses overlaid on masks
    axes[1].imshow(labeled_mask, cmap='gray', alpha=0.7)
    
    for result in ellipse_results:
        # Convert coordinates back to pixels for visualization
        center_x_pixels = result['center_x'] / pixel_size_um
        center_y_pixels = result['center_y'] / pixel_size_um
        major_axis_pixels = result['major_axis'] / pixel_size_um
        minor_axis_pixels = result['minor_axis'] / pixel_size_um
        
        # Create ellipse patch with pixel coordinates
        from matplotlib.patches import Ellipse as MPLEllipse
        ell = MPLEllipse((center_x_pixels, center_y_pixels),
                        major_axis_pixels, minor_axis_pixels,
                        angle=result['angle'],
                        fill=False, edgecolor='red', linewidth=2)
        axes[1].add_patch(ell)
        
        # Add cell ID annotation
        axes[1].text(center_x_pixels, center_y_pixels, 
                    str(result['cell_id']), 
                    color='yellow', fontsize=8, ha='center', va='center')
    
    axes[1].set_title('Fitted Ellipses')
    axes[1].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Visualization saved to: {output_path}")


def get_tiff_files(directory):
    """Get all TIFF files in a directory."""
    tiff_patterns = ['*.tif', '*.tiff', '*.TIF', '*.TIFF']
    tiff_files = []
    
    for pattern in tiff_patterns:
        tiff_files.extend(glob.glob(os.path.join(directory, pattern)))
    
    return sorted(tiff_files)


def process_single_image(image_path, output_dir=None, save_individual=True, pixel_size_um=1.0):
    """Process a single image and return ellipse results."""
    image_name = os.path.basename(image_path)
    print(f"Processing: {image_name}")
    
    # Load image
    image = load_mask_image(image_path)
    if image is None:
        return []
    
    # Preprocess mask
    labeled_mask = preprocess_mask(image)
    
    if labeled_mask.max() == 0:
        print(f"No cell instances found in {image_name}")
        return []
    
    # Fit ellipses to all cells
    ellipse_results = fit_ellipses_to_all_cells(labeled_mask, pixel_size_um)
    
    # Add image filename to results
    for result in ellipse_results:
        result['image_file'] = image_name
    
    print(f"Found {len(ellipse_results)} cells in {image_name}")
    
    # Save individual results if requested
    if save_individual and output_dir and ellipse_results:
        save_individual_results(image, labeled_mask, ellipse_results, image_name, output_dir, pixel_size_um)
    
    return ellipse_results


def save_individual_results(image, labeled_mask, ellipse_results, image_name, output_dir, pixel_size_um=1.0):
    """Save individual results for a single image."""
    # Create individual images directory
    individual_dir = os.path.join(output_dir, 'individual_images')
    os.makedirs(individual_dir, exist_ok=True)
    
    # Base name without extension
    base_name = os.path.splitext(image_name)[0]
    
    # Save individual CSV
    df_individual = pd.DataFrame(ellipse_results)
    csv_path = os.path.join(individual_dir, f"{base_name}_ellipse_parameters.csv")
    df_individual.to_csv(csv_path, index=False)
    
    # Create individual visualization
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    # Original image
    axes[0].imshow(image, cmap='gray')
    axes[0].set_title(f'Original Image: {image_name}')
    axes[0].axis('off')
    
    # Labeled masks
    axes[1].imshow(labeled_mask, cmap='tab20')
    axes[1].set_title(f'Cell Masks ({len(ellipse_results)} cells)')
    axes[1].axis('off')
    
    # Ellipses overlaid on masks
    axes[2].imshow(labeled_mask, cmap='gray', alpha=0.7)
    
    for result in ellipse_results:
        # Convert coordinates back to pixels for visualization
        center_x_pixels = result['center_x'] / pixel_size_um
        center_y_pixels = result['center_y'] / pixel_size_um
        major_axis_pixels = result['major_axis'] / pixel_size_um
        minor_axis_pixels = result['minor_axis'] / pixel_size_um
        
        # Create ellipse patch
        from matplotlib.patches import Ellipse as MPLEllipse
        ell = MPLEllipse((center_x_pixels, center_y_pixels),
                        major_axis_pixels, minor_axis_pixels,
                        angle=result['angle'],
                        fill=False, edgecolor='red', linewidth=2)
        axes[2].add_patch(ell)
        
        # Add cell ID annotation
        axes[2].text(center_x_pixels, center_y_pixels, 
                    str(result['cell_id']), 
                    color='yellow', fontsize=10, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.7))
    
    axes[2].set_title('Fitted Ellipses')
    axes[2].axis('off')
    
    plt.tight_layout()
    
    # Save visualization
    viz_path = os.path.join(individual_dir, f"{base_name}_ellipse_fit.png")
    plt.savefig(viz_path, dpi=300, bbox_inches='tight')
    plt.close()  # Close to save memory
    
    # Create individual parameter histogram
    if len(ellipse_results) > 1:  # Only create histogram if more than 1 cell
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        parameters = {
            'area': 'Area (µm²)',
            'volume': 'Volume (µm³)',
            'eccentricity': 'Eccentricity',
            'major_axis': 'Major Axis (µm)',
            'minor_axis': 'Minor Axis (µm)'
        }
        
        for i, (param, label) in enumerate(parameters.items()):
            values = [result[param] for result in ellipse_results]
            axes[i].hist(values, bins=min(10, len(values)), alpha=0.7, 
                        color='skyblue', edgecolor='black')
            axes[i].set_xlabel(label)
            axes[i].set_ylabel('Frequency')
            axes[i].set_title(f'{label} - {image_name}')
            axes[i].grid(True, alpha=0.3)
            
            # Add mean line
            mean_val = np.mean(values)
            axes[i].axvline(mean_val, color='red', linestyle='--', linewidth=2,
                           label=f'Mean: {mean_val:.2f}')
            axes[i].legend()
        
        # Hide the empty subplot
        if len(parameters) < len(axes):
            axes[-1].set_visible(False)
        
        plt.tight_layout()
        hist_path = os.path.join(individual_dir, f"{base_name}_parameter_histogram.png")
        plt.savefig(hist_path, dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Individual results saved for {image_name}")


def process_all_images(input_directory, output_dir=None, pixel_size_um=1.0):
    """Process all TIFF images in a directory."""
    tiff_files = get_tiff_files(input_directory)
    
    if not tiff_files:
        print(f"No TIFF files found in {input_directory}")
        return []
    
    print(f"Found {len(tiff_files)} TIFF files to process")
    
    all_results = []
    
    for tiff_file in tiff_files:
        results = process_single_image(tiff_file, output_dir, save_individual=True, pixel_size_um=pixel_size_um)
        all_results.extend(results)
    
    return all_results


def create_histograms(df, output_dir):
    """Create histograms of ellipse parameters."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Parameters to plot
    parameters = {
        'area': 'Area (µm²)',
        'volume': 'Volume (µm³)',
        'eccentricity': 'Eccentricity',
        'major_axis': 'Major Axis Length (µm)',
        'minor_axis': 'Minor Axis Length (µm)'
    }
    
    # Create individual histograms
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for i, (param, label) in enumerate(parameters.items()):
        axes[i].hist(df[param], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        axes[i].set_xlabel(label)
        axes[i].set_ylabel('Frequency')
        axes[i].set_title(f'Distribution of {label}')
        axes[i].grid(True, alpha=0.3)
        
        # Add statistics to the plot
        mean_val = df[param].mean()
        std_val = df[param].std()
        axes[i].axvline(mean_val, color='red', linestyle='--', linewidth=2, 
                       label=f'Mean: {mean_val:.2f}')
        axes[i].axvline(mean_val + std_val, color='orange', linestyle=':', linewidth=1,
                       label=f'+1σ: {mean_val + std_val:.2f}')
        axes[i].axvline(mean_val - std_val, color='orange', linestyle=':', linewidth=1,
                       label=f'-1σ: {mean_val - std_val:.2f}')
        axes[i].legend(fontsize=8)
    
    # Hide the empty subplot
    if len(parameters) < len(axes):
        axes[-1].set_visible(False)
    
    plt.tight_layout()
    
    # Save combined histogram
    hist_path = os.path.join(output_dir, 'ellipse_parameters_histograms.png')
    plt.savefig(hist_path, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Histograms saved to: {hist_path}")
    
    # Create summary statistics plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Box plot of all parameters (normalized)
    normalized_data = []
    labels = []
    
    for param, label in parameters.items():
        # Normalize to 0-1 scale for comparison
        data = df[param]
        normalized = (data - data.min()) / (data.max() - data.min())
        normalized_data.append(normalized)
        labels.append(param.replace('_', ' ').title())
    
    box_plot = ax.boxplot(normalized_data, labels=labels, patch_artist=True)
    
    # Color the boxes
    colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow']
    for patch, color in zip(box_plot['boxes'], colors):
        patch.set_facecolor(color)
    
    ax.set_ylabel('Normalized Values (0-1)')
    ax.set_title('Distribution Comparison of Ellipse Parameters')
    ax.grid(True, alpha=0.3)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save box plot
    box_path = os.path.join(output_dir, 'ellipse_parameters_boxplot.png')
    plt.savefig(box_path, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Box plot saved to: {box_path}")


def main():
    parser = argparse.ArgumentParser(description='Fit ellipses to cell mask instances')
    parser.add_argument('--input', '-i', 
                       default='../mitochondria_FITC/Y333 ATP6 ATP2/aligned_manual/',
                       help='Input directory containing TIFF images or single TIFF file path')
    parser.add_argument('--output', '-o', 
                       default='cell_ellipse_parameters_all.csv',
                       help='Output CSV file path')
    parser.add_argument('--output_dir', '-d',
                       default='ellipse_analysis_results',
                       help='Output directory for results and visualizations')
    parser.add_argument('--pixel_size', '-p',
                       type=float, default=1.0,
                       help='Pixel size in micrometers per pixel (default: 1.0 µm/pixel)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Check if input is a directory or file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input path does not exist: {input_path}")
        return
    
    if input_path.is_file():
        # Process single file
        print(f"Processing single file: {input_path}")
        print(f"Using pixel size: {args.pixel_size} µm/pixel")
        all_results = process_single_image(str(input_path), args.output_dir, save_individual=True, pixel_size_um=args.pixel_size)
    elif input_path.is_dir():
        # Process all files in directory
        print(f"Processing all TIFF files in directory: {input_path}")
        print(f"Using pixel size: {args.pixel_size} µm/pixel")
        all_results = process_all_images(str(input_path), args.output_dir, pixel_size_um=args.pixel_size)
    else:
        print(f"Error: Invalid input path: {input_path}")
        return
    
    if not all_results:
        print("No successful ellipse fits found")
        return
    
    # Convert results to DataFrame
    df = pd.DataFrame(all_results)
    
    # Reorder columns for better readability
    column_order = ['image_file', 'cell_id', 'center_x', 'center_y', 'major_axis', 'minor_axis', 
                   'angle', 'area', 'volume', 'eccentricity', 'contour_area']
    df = df[column_order]
    
    # Save to CSV with metadata
    output_path = os.path.join(args.output_dir, args.output)
    
    # Create metadata header
    with open(output_path, 'w') as f:
        f.write(f"# Cell ellipse analysis results\n")
        f.write(f"# Pixel size: {args.pixel_size} µm/pixel\n")
        f.write(f"# Coordinates and dimensions are in micrometers\n")
        f.write(f"# Areas are in µm², volumes are in µm³\n")
        f.write(f"# Generated on: {pd.Timestamp.now()}\n")
        f.write("#\n")
    
    # Append the dataframe
    df.to_csv(output_path, mode='a', index=False)
    print(f"Results saved to: {output_path}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"Pixel size used: {args.pixel_size} µm/pixel")
    print(f"Total images processed: {df['image_file'].nunique()}")
    print(f"Total cells processed: {len(all_results)}")
    print(f"Average cells per image: {len(all_results) / df['image_file'].nunique():.1f}")
    print(f"Average major axis: {df['major_axis'].mean():.2f} ± {df['major_axis'].std():.2f} µm")
    print(f"Average minor axis: {df['minor_axis'].mean():.2f} ± {df['minor_axis'].std():.2f} µm")
    print(f"Average eccentricity: {df['eccentricity'].mean():.3f} ± {df['eccentricity'].std():.3f}")
    print(f"Average area: {df['area'].mean():.2f} ± {df['area'].std():.2f} µm²")
    print(f"Average volume: {df['volume'].mean():.2f} ± {df['volume'].std():.2f} µm³")
    
    # Create histograms
    create_histograms(df, args.output_dir)
    
    # Print per-image summary
    print("\nPer-image summary:")
    image_summary = df.groupby('image_file').agg({
        'cell_id': 'count',
        'area': ['mean', 'std'],
        'volume': ['mean', 'std'],
        'eccentricity': ['mean', 'std'],
        'major_axis': ['mean', 'std'],
        'minor_axis': ['mean', 'std']
    }).round(2)
    
    image_summary.columns = ['cell_count', 'area_mean', 'area_std', 'volume_mean', 'volume_std',
                           'ecc_mean', 'ecc_std', 'major_mean', 'major_std', 'minor_mean', 'minor_std']
    
    summary_path = os.path.join(args.output_dir, 'per_image_summary.csv')
    image_summary.to_csv(summary_path)
    print(f"Per-image summary saved to: {summary_path}")
    print(image_summary.head())
    
    print(f"\nAnalysis complete! All results saved to: {args.output_dir}")
    print(f"Individual image results saved to: {os.path.join(args.output_dir, 'individual_images')}")
    print("\nGenerated files:")
    print(f"  - Combined results: {output_path}")
    print(f"  - Population histograms: {os.path.join(args.output_dir, 'ellipse_parameters_histograms.png')}")
    print(f"  - Box plots: {os.path.join(args.output_dir, 'ellipse_parameters_boxplot.png')}")
    print(f"  - Per-image summary: {summary_path}")
    print(f"  - Individual results: {os.path.join(args.output_dir, 'individual_images')}/")
    print(f"    └── Each image has: CSV file, visualization, and histogram")


if __name__ == "__main__":
    main()