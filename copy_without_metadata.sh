#!/bin/bash

# Script to copy files without creating ._ AppleDouble files on macOS

# Method 1: Using rsync (recommended)
copy_with_rsync() {
    local source="$1"
    local destination="$2"
    
    echo "Copying with rsync (excludes ._ files)..."
    rsync -av --exclude="._*" --exclude=".DS_Store" "$source" "$destination"
}

# Method 2: Using tar (preserves structure, no metadata)
copy_with_tar() {
    local source="$1"
    local destination="$2"
    
    echo "Copying with tar (no resource forks)..."
    cd "$(dirname "$source")"
    tar -cf - "$(basename "$source")" | (cd "$destination" && tar -xf -)
}

# Method 3: Using cp with COPYFILE_DISABLE
copy_with_cp_disabled() {
    local source="$1"
    local destination="$2"
    
    echo "Copying with cp and COPYFILE_DISABLE=1..."
    COPYFILE_DISABLE=1 cp -R "$source" "$destination"
}

# Method 4: Clean up after copying
cleanup_after_copy() {
    local destination="$1"
    
    echo "Cleaning up ._ files after copy..."
    find "$destination" -name "._*" -type f -delete
    find "$destination" -name ".DS_Store" -type f -delete
}

# Usage examples:
echo "Usage examples:"
echo "copy_with_rsync /source/path /destination/path"
echo "copy_with_tar /source/path /destination/path"
echo "copy_with_cp_disabled /source/path /destination/path"
echo "cleanup_after_copy /destination/path"

# For your specific case, you can use:
echo ""
echo "For your current directory structure, run:"
echo "cleanup_after_copy '/Volumes/ExFAT/mRNA_molecules/Y333 ATP6 ATP2/extracted_cells_conn'" 