#!/bin/bash

# Set parameters
mitograph_exe="./MitoGraph_z"
xy=0.0645
z=0.2
threshold=0.2          # 提高阈值以降低敏感性，减少碎片化
scales_min=1.0         # 最小尺度，用于检测细小结构
scales_max=3.0         # 最大尺度，用于检测较粗的管状结构
scales_count=8         # 尺度数量，更多尺度有助于更好的管状结构检测
z_adaptive=true
z_block_size=15        # 增加z_block_size以获得更平滑的结果
root_path="Y333 ATP6 ATP2/extracted_cells_z_adaptive"

echo "Script started"
echo "MitoGraph executable: $mitograph_exe"
echo "Parameters: xy=$xy, z=$z, threshold=$threshold"
echo "Scales: $scales_min to $scales_max with $scales_count levels"
echo "Z-adaptive: $z_adaptive"
echo "Z-block size: $z_block_size"
echo "Root path: $root_path"

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

# Get all subdirectories
echo "Processing $root_path"
folder_count=$(find "$root_path" -maxdepth 1 -type d | wc -l)
folder_count=$((folder_count - 1))  # Subtract 1 for the root directory itself
echo "Found $folder_count folders to process"

# Process each folder
find "$root_path" -maxdepth 1 -type d -not -path "$root_path" | while read -r folder_path; do
    folder_name=$(basename "$folder_path")
    echo "Processing $folder_path"
    
    # Build command with comprehensive parameters to improve tubular structure detection
    cmd_args=(-xy "$xy" -z "$z" -path "$folder_path" -threshold "$threshold" -scales "$scales_min" "$scales_max" "$scales_count")
    
    if [ "$z_adaptive" = true ]; then
        cmd_args+=(-z-adaptive)
        cmd_args+=(-z-block-size "$z_block_size")
    fi
    
    # Execute the command
    if "$mitograph_exe" "${cmd_args[@]}"; then
        echo "Completed processing $folder_name"
    else
        echo "Error processing $folder_path" >&2
    fi
done

echo "Script completed"