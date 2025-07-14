batch_extract_cells.py - Extracts cells based on masks, creates a **coordinate_mapping.json** for mapping.

interactive_3d_cell2_batch_fixed.py - Loads the skeleton and FISH Quant spots, the coordinates will be mapped based on coordinate_mapping.json. Calculates distance and create histograms.

## Folder Structure

```
mRNA_molecules/
├── batch_extract_cells.py          # Cell extraction script
├── interactive_3d_cell2_batch_fixed.py  # Distance analysis script
├── [Experiment Folders]/           # e.g., Y333 ATP6 ATP2/, Y333 ATP6 ATP3/, Y333 ATP6 TIM50/
│   ├── aligned_masks/              # Cell segmentation masks (input for batch_extract_cells.py)
│   ├── [gene]_spots/               # FISH spot coordinates (e.g., atp2_spots/, atp6_spots/, tim50_spots/)
│   ├── deconvolved/                # Deconvolved microscopy images
│   ├── extracted_cells/            # Output from batch_extract_cells.py
│   │   ├── coordinate_mapping.json # Coordinate mapping file (created by batch_extract_cells.py)
│   │   └── [sample_directories]/   # Individual cell data (e.g., yWL333_cy3_ATP2_cy5_ATP6MS2_1_s1/)
│   │       ├── *.vtk              # 3D cell structure files
│   │       ├── *.coo              # Coordinate files
│   │       ├── *.gnet             # Network/skeleton files
│   │       └── *.pts              # Point cloud files
│   └── interactive_batch_results/  # Output from interactive_3d_cell2_batch_fixed.py
│       └── *.csv                  # Distance analysis results and histograms
```

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