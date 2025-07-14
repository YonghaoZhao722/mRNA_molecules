batch_extract_cells.py - Extracts cells based on masks, creates a **coordinate_mapping.json** for mapping.

interactive_3d_cell2_batch_fixed.py - Loads the skeleton and FISH Quant spots, the coordinates will be mapped based on coordinate_mapping.json. Calculates distance and create histograms.

## File Placement Requirements

### Before Running batch_extract_cells.py:
```
[Experiment Folder]/               # e.g., Y333 ATP6 ATP2/
├── aligned_masks/                 # REQUIRED: Place your cell segmentation mask files here
│   ├── mask_file_1.tif           # Cell segmentation masks
│   ├── mask_file_2.tif
│   └── ...
└── deconvolved/                   # OPTIONAL: Deconvolved images (if needed)
    ├── image_file_1.tif
    └── ...
```

### Before Running interactive_3d_cell2_batch_fixed.py:
```
[Experiment Folder]/               # e.g., Y333 ATP6 ATP2/
├── extracted_cells/               # OUTPUT from batch_extract_cells.py (auto-generated)
│   ├── coordinate_mapping.json   # Mapping file (auto-generated)
│   └── [sample_dirs]/            # Cell data (auto-generated)
├── atp2_spots/                   # REQUIRED: Place your FISH spot coordinate files here
│   ├── spot_file_1.txt
│   ├── spot_file_2.txt
│   └── ...
├── atp6_spots/                   # REQUIRED: Place your second gene spot files here
│   ├── spot_file_1.txt
│   └── ...
└── (other gene)_spots/           # REQUIRED: Additional gene spot directories as needed
    └── ...
```