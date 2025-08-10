#!/bin/bash

# Set parameters
mitograph_exe="./MitoGraph_no_flip"
xy=0.0645
z=0.2
threshold=0.3
scales_min=1.1
scales_max=1.5
scales_count=8
z_adaptive=true
z_block_size=15
root_path="Y333 ATP6 ATP2/extracted_cells_dw_30"
csv_file="Y333 ATP6 ATP2/interactive_batch_results/atp6_filtered_spots_exceed_threshold_0.5_extracted_cells_dw_30.csv"
temp_base_dir="temp_mito_processing"

echo "Script started"
echo "MitoGraph executable: $mitograph_exe"
echo "Parameters: xy=$xy, z=$z, threshold=$threshold"
echo "Scales: $scales_min to $scales_max with $scales_count levels"
echo "Z-adaptive: $z_adaptive"
echo "Z-block size: $z_block_size"
echo "Root path: $root_path"
echo "CSV filter file: $csv_file"
echo "Temporary processing directory: $temp_base_dir"

# Check if MitoGraph executable exists
if [ ! -e "$mitograph_exe" ]; then
    echo "Error: MitoGraph executable not found at $mitograph_exe" >&2
    exit 1
fi

# Check if root path exists
if [ ! -d "$root_path" ]; then
    echo "Error: Root path not found: $root_path" >&2
    exit 1
fi

# Check if CSV file exists
if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file not found: $csv_file" >&2
    exit 1
fi

# Function to convert Cell_X to cell_XXX format
convert_cell_name() {
    local cell_name="$1"
    local cell_number=$(echo "$cell_name" | sed 's/Cell_//')
    printf "cell_%03d" "$cell_number"
}

# Function to copy only tif files to temporary directory
copy_tif_files() {
    local source_dir="$1"
    local temp_dir="$2"
    local cell_prefix="$3"
    
    # Copy only .tif files for this cell (exclude macOS resource fork files)
    find "$source_dir" -name "${cell_prefix}.tif" -not -name "._*" -exec cp {} "$temp_dir/" \;
    
    # Also copy configuration files if they exist
    if [ -f "$source_dir/mitograph.config" ]; then
        cp "$source_dir/mitograph.config" "$temp_dir/"
    fi
    if [ -f "$source_dir/coordinate_mapping.json" ]; then
        cp "$source_dir/coordinate_mapping.json" "$temp_dir/"
    fi
}

# Function to copy all generated results back to original directory
copy_results_back() {
    local temp_dir="$1"
    local target_dir="$2"
    local cell_prefix="$3"
    
    # Copy back all generated files (exclude the original .tif input files)
    # Suppress extended attribute warnings by redirecting stderr
    find "$temp_dir" -name "${cell_prefix}*" \( -name "*.gnet" -o -name "*.mitograph" -o -name "*.coo" -o -name "*.vtk" -o -name "*.png" -o -name "*.txt" \) -exec cp {} "$target_dir/" \; 2>/dev/null || true
}

# Create main temporary directory
if [ -d "$temp_base_dir" ]; then
    echo "Removing existing temporary directory: $temp_base_dir"
    rm -rf "$temp_base_dir"
fi
mkdir -p "$temp_base_dir"

# Extract unique images from CSV file
echo "Extracting images and cells from CSV file..."

unique_images=$(tail -n +2 "$csv_file" | cut -d',' -f1 | sort | uniq)
image_count=$(echo "$unique_images" | wc -l)
echo "Found $image_count unique images to process"

# Process each image and its cells
processed_images=0
skipped_images=0
total_processed_cells=0

# Convert unique_images to array for proper variable scope
IFS=$'\n' read -d '' -r -a images_array <<< "$unique_images"

for image_name in "${images_array[@]}"; do
    # Get all cells for this image
    cells_for_image=$(tail -n +2 "$csv_file" | awk -F',' -v img="$image_name" '$1==img {print $2}' | sort | uniq | tr '\n' ' ')
    source_folder="$root_path/$image_name"
    
    if [ -d "$source_folder" ]; then
        echo "Processing image: $image_name"
        echo "  Cells to process: $cells_for_image"
        
        # Create temporary directory for this image
        temp_dir="$temp_base_dir/$image_name"
        mkdir -p "$temp_dir"
        
        # Copy required cell tif files and configuration files
        cells_copied=0
        for cell_name in $cells_for_image; do
            cell_prefix=$(convert_cell_name "$cell_name")
            
            # Check if cell tif file exists before copying
            if [ -f "$source_folder/${cell_prefix}.tif" ]; then
                # Copy only the tif file for this cell
                copy_tif_files "$source_folder" "$temp_dir" "$cell_prefix"
                echo "    Copied tif file for $cell_name ($cell_prefix)"
                cells_copied=$((cells_copied + 1))
            else
                echo "    Warning: No tif file found for $cell_name ($cell_prefix)"
            fi
        done
        
        # Clean up any macOS resource fork files that might have been copied
        find "$temp_dir" -name "._*" -delete 2>/dev/null || true
        
        if [ "$cells_copied" -eq 0 ]; then
            echo "  No valid cells found for $image_name, skipping"
            skipped_images=$((skipped_images + 1))
            rm -rf "$temp_dir"
            continue
        fi
        
        echo "  Successfully copied files for $cells_copied cells"
        
        # Build command with comprehensive parameters
        cmd_args=(-xy "$xy" -z "$z" -path "$temp_dir" -threshold "$threshold" -scales "$scales_min" "$scales_max" "$scales_count")
        
        if [ "$z_adaptive" = true ]; then
            cmd_args+=(-z-adaptive)
            cmd_args+=(-z-block-size "$z_block_size")
        fi
        
        # Execute MitoGraph on temporary directory
        echo "  Running MitoGraph on $temp_dir..."
        if "$mitograph_exe" "${cmd_args[@]}"; then
            echo "  MitoGraph processing completed for $image_name"
            
            # Copy results back to original directory for each cell
            for cell_name in $cells_for_image; do
                cell_prefix=$(convert_cell_name "$cell_name")
                copy_results_back "$temp_dir" "$source_folder" "$cell_prefix"
            done
            echo "  Results copied back to $source_folder"
            
            processed_images=$((processed_images + 1))
            total_processed_cells=$((total_processed_cells + cells_copied))
        else
            echo "  Error: MitoGraph processing failed for $temp_dir" >&2
        fi
        
        # Clean up temporary directory
        rm -rf "$temp_dir"
        
    else
        echo "Warning: Source folder not found: $source_folder"
        skipped_images=$((skipped_images + 1))
    fi
done

# Clean up main temporary directory
if [ -d "$temp_base_dir" ]; then
    rm -rf "$temp_base_dir"
    echo "Cleaned up temporary directory: $temp_base_dir"
fi

echo "Script completed"
echo "Total processed images: $processed_images"
echo "Total skipped images: $skipped_images"
echo "Total processed cells: $total_processed_cells"
