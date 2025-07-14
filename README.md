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

### Workflow
1. **batch_extract_cells.py**: Processes aligned_masks/ to extract individual cells into extracted_cells/ and generates coordinate_mapping.json
2. **interactive_3d_cell2_batch_fixed.py**: Uses extracted_cells/ data, [gene]_spots/, and coordinate_mapping.json to perform distance analysis and save results to interactive_batch_results/