import os
import numpy as np
from bigfish import stack, detection
from pathlib import Path
import tifffile
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def visualize_spots(img, spots, save_path, threshold):
    """
    Visualize detected spots on maximum intensity projection and individual z-slices
    
    Args:
        img (np.ndarray): 3D image array
        spots (np.ndarray): Detected spot coordinates
        save_path (str): Path to save visualization
        threshold (float): Threshold used for spot detection
    """
    # Create figure with subplots
    n_z = min(5, img.shape[0])  # Show up to 5 z-slices
    fig = plt.figure(figsize=(15, 5 + 3 * ((n_z + 1) // 2)))
    
    # Plot maximum intensity projection
    plt.subplot(((n_z + 1) // 2 + 1), 2, 1)
    mip = np.max(img, axis=0)
    plt.imshow(mip, cmap='gray')
    
    # Plot spots on MIP
    if len(spots) > 0:
        # Get unique z positions
        unique_z = np.unique(spots[:, 0])
        # Create colormap for z positions
        colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_z)))
        z_color_dict = dict(zip(unique_z, colors))
        
        for spot in spots:
            z, y, x = spot
            circle = Circle((x, y), radius=3, fill=False, 
                          color=z_color_dict[z], linewidth=1)
            plt.gca().add_patch(circle)
    
    plt.title(f'Maximum Intensity Projection\n{len(spots)} spots detected (threshold={threshold:.3f})')
    plt.colorbar()
    
    # Plot individual z-slices
    z_step = max(1, img.shape[0] // n_z)
    for i, z in enumerate(range(0, img.shape[0], z_step)):
        if i >= n_z:
            break
        plt.subplot(((n_z + 1) // 2 + 1), 2, i + 3)
        plt.imshow(img[z], cmap='gray')
        
        # Plot spots in this z-plane
        if len(spots) > 0:
            z_spots = spots[spots[:, 0] == z]
            for spot in z_spots:
                _, y, x = spot
                circle = Circle((x, y), radius=3, fill=False, 
                              color='red', linewidth=1)
                plt.gca().add_patch(circle)
        
        plt.title(f'Z-slice {z}')
        plt.colorbar()
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def process_mrna_spots(input_path, voxel_size=(100, 100, 100), spot_radius=(200, 200, 200)):
    """
    Process 3D mRNA FISH data to detect spots and extract their coordinates
    
    Args:
        input_path (str): Path to the directory containing FISH images
        voxel_size (tuple): Size of a voxel in nanometers (z, y, x)
        spot_radius (tuple): Expected radius of spots in nanometers (z, y, x)
    """
    # Get all .tif files in the directory
    image_files = list(Path(input_path).glob('*.tif'))
    
    # Filter out DIC and DAPI images
    image_files = [f for f in image_files if 'DIC' not in f.name and 'DAPI' not in f.name]
    print(f"Found {len(image_files)} mRNA FISH images to process")
    
    results = []
    
    for img_path in image_files:
        print(f"\nProcessing {img_path.name}")
        
        # Read the image
        img = tifffile.imread(str(img_path))
        
        # Ensure the image is in the correct format (3D)
        if len(img.shape) == 2:
            img = np.expand_dims(img, axis=0)
        
        # Convert to float32 and normalize if needed
        img = img.astype(np.float32)
        if img.max() > 1:
            img = img / img.max()
        # Detect spots using LoG filter with automatic threshold detection
        spots, threshold = detection.detect_spots(
            images=img,
            remove_duplicate=True,
            return_threshold=True,
            voxel_size=voxel_size,
            spot_radius=spot_radius
        )
        print(f"Detected spots using threshold: {threshold}")
        print(f"Number of spots detected: {len(spots)}")
        # Visualize results
        vis_path = str(img_path).rsplit('.', 1)[0] + '_spots.png'
        visualize_spots(img, spots, vis_path, threshold)
        print(f"Visualization saved to: {vis_path}")
        
        # Extract coordinates
        if len(spots) > 0:
            print(f"Found {len(spots)} spots")
            for spot in spots:
                results.append({
                    'image_name': img_path.name,
                    'z': spot[0],
                    'y': spot[1],
                    'x': spot[2]
                })
        else:
            print("No spots detected in this image")
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)
    
    # Save results
    output_file = os.path.join(input_path, 'spot_coordinates.csv')
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    return df

if __name__ == "__main__":
    # Process both directories
    paths = [
        r"F:\atp\Y333 ATP6 ATP2",
        r"F:\atp\Y333 ATP6 ATP3"
    ]
    
    for path in paths:
        print(f"\nProcessing directory: {path}")
        process_mrna_spots(path) 