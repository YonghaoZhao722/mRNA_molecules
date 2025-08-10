#!/bin/bash

# Configuration
PSF_FILE='/Volumes/ExFAT/mitochondria_FITC/generated_psf.tif'
INPUT_DIR='/Volumes/ExFAT/mitochondria_FITC/Y333 ATP6 ATP3/FITC'
OUTPUT_DIR='/Volumes/ExFAT/mitochondria_FITC/Y333 ATP6 ATP3/dw_30'

# --- Script Logic ---
mkdir -p "$OUTPUT_DIR"
echo "Starting batch processing..."
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "PSF File: $PSF_FILE"

# Iterate over all .TIF files in the input directory
# Note the change from *.tif to *.TIF
for input_file in "$INPUT_DIR"/*.TIF; do
  if [ -f "$input_file" ]; then
    filename=$(basename -- "$input_file")
    output_file="$OUTPUT_DIR/dw_$filename"

    echo "-------------------------------------"
    echo "Processing file: $filename"
    echo "Full input path: $input_file"
    echo "Full output path: $output_file"
    
    dw --iter 30 --out "$output_file" "$input_file" "$PSF_FILE"
    
  else
    # This block will now only execute if no .TIF files are found
    echo "No .TIF files found in $INPUT_DIR"
    # To prevent printing this multiple times, you could add an 'exit 1' or move this message outside the loop.
  fi
done

echo "-------------------------------------"
echo "Batch processing complete."