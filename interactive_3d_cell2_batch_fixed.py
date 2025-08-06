import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.optimize import minimize
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import threading
import tkinter.simpledialog as simpledialog
import json
import tifffile
from skimage.measure import find_contours
import vtk

class Interactive3DBatchVisualizerFixed:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive 3D Cell Batch Visualization (Fixed)")
        self.root.geometry("1500x950")

        # Path settings
        self.base_dir = r'Y333 ATP6 ATP2'
        self.skeleton_type = 'extracted_cells'
        self.skeleton_root = os.path.join(self.base_dir, self.skeleton_type)
        self.channel = 'atp6_filtered'
        self.spots_root = os.path.join(self.base_dir, f'{self.channel}_spots')
        self.mask_root = os.path.join(self.base_dir, 'aligned_masks')

        # Data structures
        self.available_images = {}  # {image_name: {cell_name: {spots_data, skeleton_file}}}
        self.coordinate_mappings = {} # To store coordinate mappings from json
        self.selected_image = tk.StringVar()
        self.selected_cell = tk.StringVar()
        self.current_cell = None
        self.current_image = None

        # Data variables
        self.original_spots = None
        self.original_skeleton = None
        self.current_skeleton = None
        self.current_spots = None
        self.current_outline = None  # New: store cell outline data
        self.distances = None
        self.current_skeleton_translation = np.array([0, 0, 0])
        
        # VTK polylines for proper skeleton visualization
        self.skeleton_segments = None  # Store polylines from VTK file
        
        # Control variables
        self.use_y_flip = tk.BooleanVar(value=False)
        self.use_z_flip = tk.BooleanVar(value=False)
        self.auto_compare_yflip = tk.BooleanVar(value=True)
        self.auto_translate_skeleton = tk.BooleanVar(value=False)
        self.show_cell_outline = tk.BooleanVar(value=False)  # New: control cell outline visibility
        self.skeleton_as_lines = tk.BooleanVar(value=True)  # New: control skeleton display mode (lines vs points)

        # Create interface
        self.create_interface()
        # Scan all images and cells
        self.scan_available_images_and_cells()
        # Auto load first image and cell
        if self.available_images:
            first_image = list(self.available_images.keys())[0]
            self.selected_image.set(first_image)
            self.on_image_changed()

    def load_spots_fishquant_analyze_method(self, file_path, cell_number=1, flip_y=True, mapping_data=None, silent=False):
        """
        Complete copy of the correct load_spots_fishquant method from analyze_alignment.py
        Fix: Properly handle cases where some cells have no spots data
        """
        # Read file to find SPOTS section
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find pixel size information
        pixel_xy, pixel_z = None, None
        for i, line in enumerate(lines):
            if line.startswith('Pix-XY'):
                # Pixel size values are on the next line
                if i + 1 < len(lines):
                    next_line = lines[i + 1].strip()
                    parts = next_line.split('\t')
                    if len(parts) >= 2:
                        pixel_xy = float(parts[0])
                        pixel_z = float(parts[1])
                break
        
        if not silent:
            print(f"FISH-QUANT pixel size: XY={pixel_xy}nm, Z={pixel_z}nm")
        
        # Find SPOTS data for the specified cell
        target_cell = f"CELL\tCell_{cell_number}"
        cell_found = False
        spots_start = -1
        spots_end = -1
        
        for i, line in enumerate(lines):
            line_stripped = line.strip()
            
            # Found target cell marker
            if line_stripped == target_cell:
                cell_found = True
                if not silent:
                    print(f"Found target cell: Cell_{cell_number}")
                continue
            
            # If target cell found and SPOTS marker encountered (indicating this cell has spots data)
            if cell_found and line.startswith('Pos_Y'):
                spots_start = i
                continue
            
            # If target cell and spots start position found, next CELL marker indicates end
            if cell_found and spots_start != -1 and line.startswith('CELL'):
                spots_end = i
                break
                
            # Critical fix: If target cell found but spots start position not found,
            # encountering next CELL marker means target cell has no spots data
            if cell_found and spots_start == -1 and line.startswith('CELL'):
                if not silent:
                    print(f"Cell_{cell_number} has no spots data (next cell encountered before finding spots start position)")
                raise ValueError(f"Cell_{cell_number} has no spots data")
        
        # If no next CELL marker found, this is the last cell
        if cell_found and spots_start != -1 and spots_end == -1:
            spots_end = len(lines)
        
        if not cell_found:
            raise ValueError(f"Cell_{cell_number} data not found")
        
        if spots_start == -1:
            raise ValueError(f"Cell_{cell_number} has no spots data")
        
        # Read spots data for this cell
        spots_data = []
        for i in range(spots_start + 1, spots_end):
            line = lines[i].strip()
            if line and not line.startswith('X_POS') and not line.startswith('Y_POS'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        y, x, z = float(parts[0]), float(parts[1]), float(parts[2])
                        spots_data.append([y, x, z])
                    except ValueError:
                        continue
        
        coords = np.array(spots_data)
        
        # If Y-axis flip is needed
        if flip_y and mapping_data:
            # Use crop region information from mapping data for precise flip
            cell_name = f"cell_{cell_number:03d}"
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                y_start = crop_info['y_start']  # Pixel coordinates
                y_end = crop_info['y_end']      # Pixel coordinates
                
                # Convert to nanometer coordinates (pixel * pixel size)
                y_start_nm = y_start * pixel_xy
                y_end_nm = y_end * pixel_xy
                
                # Flip based on crop region: new_y = (y_start + y_end) - old_y
                flip_center = y_start_nm + y_end_nm
                coords[:, 0] = flip_center - coords[:, 0]
                
                if not silent:
                    print(f"Y-axis flip: Based on crop region y_start={y_start_nm:.1f}nm, y_end={y_end_nm:.1f}nm")
            else:
                if not silent:
                    print(f"Warning: Mapping information for {cell_name} not found, skipping Y-axis flip")
        elif flip_y:
            if not silent:
                print("Warning: Mapping data needed for precise Y-axis flip, using simple flip")
            # Simple flip: Based on current cell spots range
            if len(coords) > 0:
                y_min, y_max = coords[:, 0].min(), coords[:, 0].max()
                flip_center = y_min + y_max
                coords[:, 0] = flip_center - coords[:, 0]
                if not silent:
                    print(f"Y-axis flip: Based on cell range, flip center={flip_center/2:.1f}nm")
        
        if not silent:
            print(f"Cell_{cell_number} spots data loading completed: {len(coords)} points")
            print(f"Coordinate range: X=[{coords[:,1].min():.1f}, {coords[:,1].max():.1f}], Y=[{coords[:,0].min():.1f}, {coords[:,0].max():.1f}], Z=[{coords[:,2].min():.1f}, {coords[:,2].max():.3f}]")
        
        return coords, pixel_xy, pixel_z

    def load_skeleton_txt_analyze_method(self, file_path, mapping_data=None, cell_name="cell_001", pixel_size_xy=0.0645, silent=False):
        """
        Complete copy of the correct load_skeleton_txt method from analyze_alignment.py
        """
        df = pd.read_csv(file_path, sep='\t')
        # Extract coordinate columns (x, y, z) - these are coordinates relative to the cropped region
        coords = df[['x', 'y', 'z']].values
        if not silent:
            print(f"Skeleton data loading completed: {len(coords)} points")
            print(f"Relative coordinate range: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
        
        # If mapping data is provided, convert to absolute coordinates
        if mapping_data and cell_name:
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                x_offset = crop_info['x_offset']
                y_offset = crop_info['y_offset']
                
                # Convert pixel offset to micrometers and add to relative coordinates
                offset_x_um = x_offset * pixel_size_xy
                offset_y_um = y_offset * pixel_size_xy
                
                coords[:, 0] += offset_x_um  # Add offset to X coordinate
                coords[:, 1] += offset_y_um  # Add offset to Y coordinate
                # Z coordinate doesn't need offset as it's the depth direction of 3D image
                
                if not silent:
                    print(f"Applied coordinate offset: X+{offset_x_um:.3f}μm (pixel {x_offset}), Y+{offset_y_um:.3f}μm (pixel {y_offset})")
                    print(f"Absolute coordinate range: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
            else:
                if not silent:
                    print(f"Warning: {cell_name} not found in mapping file")
        else:
            if not silent:
                if mapping_data:
                    print(f"Warning: mapping file does not exist")
                print("Using relative coordinates (no offset applied)")
        
        return coords 

    def scan_available_images_and_cells(self):
        """Scan all images and cells, match spot files with image folders through index and field-of-view double matching
        Fix: Accurately detect which cells actually have spots data"""
        self.available_images = {}
        if not os.path.exists(self.skeleton_root):
            print(f"Skeleton root not found: {self.skeleton_root}")
            return
            
        # Get all spot files
        spot_files = [f for f in os.listdir(self.spots_root) if f.endswith('.txt') and '_spots' in f]
        
        for image_folder in os.listdir(self.skeleton_root):
            image_path = os.path.join(self.skeleton_root, image_folder)
            if not os.path.isdir(image_path):
                continue

            # Load coordinate mapping
            mapping_file = os.path.join(image_path, 'coordinate_mapping.json')
            if os.path.exists(mapping_file):
                with open(mapping_file, 'r') as f:
                    self.coordinate_mappings[image_folder] = json.load(f)
            else:
                print(f"Warning: coordinate_mapping.json not found in {image_path}")
                continue

            # Parse index and field-of-view
            parts = image_folder.split('_')
            if len(parts) < 2:
                continue
            index = parts[-2]
            fov = parts[-1]
            
            # Match spot file
            matched_spot = None
            for spot_file in spot_files:
                if f'_{index}_' in spot_file and f'_{fov}_' in spot_file and '_spots' in spot_file:
                    matched_spot = os.path.join(self.spots_root, spot_file)
                    break
            
            if not matched_spot:
                print(f"Spot file not found for {image_folder} (index={index}, fov={fov})")
                continue

            # Scan skeleton files
            skeleton_files = [f for f in os.listdir(image_path) if f.startswith('cell_') and f.endswith('.txt')]
            image_cells = {}
            
            for skel_file in skeleton_files:
                try:
                    cell_num = int(skel_file.split('_')[1].split('.')[0])
                    cell_name = f"Cell_{cell_num}"
                    
                    # Strictly check if this cell actually has spots data in the spot file
                    try:
                        # Use the fixed function to strictly verify spots data
                        test_coords, _, _ = self.load_spots_fishquant_analyze_method(
                            matched_spot, 
                            cell_number=cell_num, 
                            flip_y=False,  # Don't flip during testing
                            mapping_data=None,  # Don't use mapping during testing
                            silent=True  # Silent mode
                        )
                        
                        # If successfully loaded with data, add to available list
                        if len(test_coords) > 0:
                            image_cells[cell_name] = {
                                'spots_file': matched_spot,
                                'skeleton_file': os.path.join(image_path, skel_file)
                            }
                            print(f"Found {image_folder}-{cell_name} (verified with {len(test_coords)} spots)")
                        else:
                            print(f"Skipping {image_folder}-{cell_name}: spots data is empty")
                            
                    except ValueError as e:
                        # If "no spots data" exception is thrown, skip this cell
                        print(f"Skipping {image_folder}-{cell_name}: {str(e)}")
                    except Exception as e:
                        print(f"Error checking spots file for {cell_name}: {e}")
                        
                except Exception as e:
                    print(f"Failed to parse skeleton filename {skel_file}: {e}")
            
            if image_cells:
                self.available_images[image_folder] = image_cells
        
        # Update image dropdown
        if hasattr(self, 'image_combo'):
            sorted_image_names = sorted(list(self.available_images.keys()))
            self.image_combo['values'] = sorted_image_names
            if self.available_images:
                self.selected_image.set(sorted_image_names[0])

    def find_mask_file_for_image(self, image_name):
        """Find the corresponding mask file for a given image name"""
        if not os.path.exists(self.mask_root):
            return None
        
        # Parse image name to get prefix and field index
        parts = image_name.split('_')
        if len(parts) < 2:
            return None
        
        # Extract index and field-of-view from image name
        index = parts[-2]  # e.g., "1", "2", "3"
        fov = parts[-1]    # e.g., "s1", "s2", "s3"
        
        # Reconstruct the prefix (everything before the last two parts)
        prefix = '_'.join(parts[:-2])
        
        # Look for mask file with pattern: prefix_DIC_index_fov.tif
        mask_filename = f"{prefix}_DIC_{index}_{fov}.tif"
        mask_path = os.path.join(self.mask_root, mask_filename)
        
        if os.path.exists(mask_path):
            return mask_path
        
        # Try alternative pattern: prefix_index_DIC_fov.tif
        mask_filename = f"{prefix}_{index}_DIC_{fov}.tif"
        mask_path = os.path.join(self.mask_root, mask_filename)
        
        if os.path.exists(mask_path):
            return mask_path
        
        return None

    def load_cell_outline_from_mask(self, image_name, cell_name, pixel_xy):
        """Load cell outline from mask file and convert to proper coordinates"""
        try:
            # Find the mask file for this image
            mask_file = self.find_mask_file_for_image(image_name)
            if not mask_file:
                print(f"Warning: No mask file found for {image_name}")
                return None
            
            # Load mask data
            mask = tifffile.imread(mask_file)
            
            # Get cell information from mapping data
            cell_number = int(cell_name.split("_")[1])
            cell_id_str = f"{cell_number:03d}"
            mapping_cell_name = f'cell_{cell_id_str}'
            mapping_data = self.coordinate_mappings.get(image_name, {})
            
            if mapping_cell_name not in mapping_data:
                print(f"Warning: No mapping data found for {mapping_cell_name}")
                return None
            
            cell_mapping = mapping_data[mapping_cell_name]
            original_mask_label = cell_mapping.get('original_mask_label', cell_number)
            
            # Extract the specific cell mask
            cell_mask = (mask == original_mask_label).astype(np.uint8)
            
            # Find contours of the cell
            contours = find_contours(cell_mask, level=0.5)
            
            if len(contours) == 0:
                print(f"Warning: No contours found for {mapping_cell_name}")
                return None
            
            # Use the largest contour (main cell boundary)
            main_contour = max(contours, key=len)
            
            # Convert contour coordinates to micrometers
            # Contour coordinates are in (y, x) format from find_contours
            # NOTE: Outline is extracted from original mask, so it's already in absolute coordinates
            # unlike skeleton which is in relative coordinates and needs offset
            
            pixel_size_um = pixel_xy / 1000.0  # Convert nm to μm
            
            # Convert contour to (x, y) coordinates in micrometers (absolute coordinates)
            outline_coords = np.zeros((len(main_contour), 2))
            outline_coords[:, 0] = main_contour[:, 1] * pixel_size_um  # x coordinates
            outline_coords[:, 1] = main_contour[:, 0] * pixel_size_um  # y coordinates
            
            # DO NOT apply coordinate offset - outline is already in absolute coordinates!
            # The offset is only for skeleton data which is in relative coordinates
            
            # Apply Y-flip if enabled (same logic as spots data)
            if self.use_y_flip.get() and mapping_cell_name in mapping_data:
                crop_info = mapping_data[mapping_cell_name]['crop_region']
                y_start = crop_info['y_start']  # Pixel coordinates
                y_end = crop_info['y_end']      # Pixel coordinates
                
                # Convert to nanometer coordinates (pixel * pixel size)
                y_start_nm = y_start * pixel_xy
                y_end_nm = y_end * pixel_xy
                
                # Flip based on crop region: new_y = (y_start + y_end) - old_y
                flip_center_nm = y_start_nm + y_end_nm
                flip_center_um = flip_center_nm / 1000.0  # Convert to μm
                outline_coords[:, 1] = flip_center_um - outline_coords[:, 1]  # Both in μm units
            
            print(f"Loaded cell outline: {len(outline_coords)} points")
            print(f"Outline range: X=[{outline_coords[:,0].min():.3f}, {outline_coords[:,0].max():.3f}], Y=[{outline_coords[:,1].min():.3f}, {outline_coords[:,1].max():.3f}]")
            
            return outline_coords
            
        except Exception as e:
            print(f"Error loading cell outline: {e}")
            return None

    def get_outline_z_coordinate(self):
        """Get the Z coordinate for the cell outline (center of the cell in Z)"""
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            z_center = np.mean(self.current_skeleton[:, 2])  # use skeleton z center
        elif self.current_spots is not None and len(self.current_spots) > 0:
            z_center = np.mean(self.current_spots[:, 2])
        else:
            z_center = 0.0
        return z_center

    def load_cell_data_analyze_method(self, image_name, cell_name):
        """
        Load cell data using the correct method verified by analyze_alignment.py
        """
        if image_name not in self.available_images or cell_name not in self.available_images[image_name]:
            print(f"Cell {cell_name} in image {image_name} not available")
            return
        
        try:
            self.current_image = image_name
            self.current_cell = cell_name
            cell_info = self.available_images[image_name][cell_name]
            
            # Get cell number and mapping data
            cell_number = int(cell_name.split("_")[1])
            cell_id_str = f"{cell_number:03d}"
            mapping_cell_name = f'cell_{cell_id_str}'
            mapping_data = self.coordinate_mappings.get(image_name, {})
            
            print(f"=== Loading {image_name} - {cell_name} ===")
            
            # 1. Load spots data (using analyze_alignment.py method)
            spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                cell_info['spots_file'], 
                cell_number=cell_number, 
                flip_y=self.use_y_flip.get(), 
                mapping_data=mapping_data
            )
            
            # Store pixel size for outline loading
            self.current_pixel_xy = pixel_xy
            
            # 2. Load skeleton data (using analyze_alignment.py method)
            skeleton_coords = self.load_skeleton_txt_analyze_method(
                cell_info['skeleton_file'], 
                mapping_data=mapping_data,
                cell_name=mapping_cell_name,
                pixel_size_xy=pixel_xy/1000  # Convert to micrometers
            )
            
            # 2.1 Load skeleton polylines from VTK file
            skeleton_dir = os.path.dirname(cell_info['skeleton_file'])
            vtk_polylines = self.load_skeleton_from_vtk(skeleton_dir, cell_name)
            
            # 3. Load cell outline data
            outline_coords = self.load_cell_outline_from_mask(image_name, cell_name, pixel_xy)
            
            # 4. Coordinate transformation (following analyze_alignment.py logic exactly)
            # Rearrange to (x, y, z) format to match skeleton
            spots_nm_xyz = spots_coords[:, [1, 0, 2]]  # from (y, x, z) to (x, y, z)
            
            # Convert to micrometer units
            spots_um_xyz = spots_nm_xyz / 1000.0
            
            print(f"Converted spots coordinate range: X=[{spots_um_xyz[:,0].min():.3f}, {spots_um_xyz[:,0].max():.3f}], Y=[{spots_um_xyz[:,1].min():.3f}, {spots_um_xyz[:,1].max():.3f}], Z=[{spots_um_xyz[:,2].min():.3f}, {spots_um_xyz[:,2].max():.3f}]")
            
            skeleton_center = np.mean(skeleton_coords, axis=0)
            spots_center = np.mean(spots_um_xyz, axis=0)
            center_diff = skeleton_center - spots_center
            center_distance = np.linalg.norm(center_diff)
            
            print(f"Center correction analysis:")
            print(f"  Skeleton center: X={skeleton_center[0]:.3f}, Y={skeleton_center[1]:.3f}, Z={skeleton_center[2]:.3f}")
            print(f"  Spots center: X={spots_center[0]:.3f}, Y={spots_center[1]:.3f}, Z={spots_center[2]:.3f}")
            print(f"  Center distance: {center_distance:.3f}μm")
            
            # Apply center correction
            spots_corrected = spots_um_xyz - center_diff
            
            # Save two versions of data to match analyze_alignment.py behavior
            self.spots_before_correction = spots_um_xyz.copy()  # Before center correction (for chart display, matching analyze_alignment.py)
            self.spots_after_correction = spots_corrected.copy()  # After center correction (for final statistics)
            
            # Use pre-correction data by default for visualization (consistent with analyze_alignment.py charts)
            self.current_spots = self.spots_before_correction.copy()
            self.original_spots = self.spots_before_correction.copy()
            self.current_skeleton = skeleton_coords.copy()
            self.original_skeleton = skeleton_coords.copy()
            
            # Apply coordinate transformation to VTK polylines
            self.skeleton_segments = self.apply_coordinate_transform_to_polylines(
                vtk_polylines,
                mapping_data=mapping_data,
                cell_name=mapping_cell_name,
                pixel_size_xy=pixel_xy/1000
            )
            
            # Process outline data for 3D visualization
            if outline_coords is not None:
                # Create 3D outline coordinates using the center Z of the cell
                z_center = self.get_outline_z_coordinate()
                outline_3d = np.zeros((len(outline_coords), 3))
                outline_3d[:, 0] = outline_coords[:, 0]  # X coordinates
                outline_3d[:, 1] = outline_coords[:, 1]  # Y coordinates  
                outline_3d[:, 2] = z_center               # Z coordinate (center level)
                
                # Close the outline loop for better visualization
                if len(outline_3d) > 0:
                    outline_3d = np.vstack([outline_3d, outline_3d[0]])  # Add first point at the end to close the loop
                
                self.current_outline = outline_3d
                print(f"Cell outline prepared for 3D visualization at Z={z_center:.3f}μm")
            else:
                self.current_outline = None
            
            print(f"Spots range before correction (for visualization): X=[{self.current_spots[:,0].min():.3f}, {self.current_spots[:,0].max():.3f}], Y=[{self.current_spots[:,1].min():.3f}, {self.current_spots[:,1].max():.3f}], Z=[{self.current_spots[:,2].min():.3f}, {self.current_spots[:,2].max():.3f}]")
            print(f"Spots range after correction: X=[{spots_corrected[:,0].min():.3f}, {spots_corrected[:,0].max():.3f}], Y=[{spots_corrected[:,1].min():.3f}, {spots_corrected[:,1].max():.3f}], Z=[{spots_corrected[:,2].min():.3f}, {spots_corrected[:,2].max():.3f}]")
            
            # Update window title
            self.root.title(f"Interactive 3D Cell Batch Visualization (Fixed) - {image_name} - {cell_name}")
            
            print(f"Successfully loaded {image_name}-{cell_name} data:")
            print(f"  Spots count: {len(self.original_spots)}")
            print(f"  Skeleton points count: {len(self.original_skeleton)}")
            
            
        except Exception as e:
            print(f"Failed to load {image_name}-{cell_name} data: {e}")
            import traceback
            traceback.print_exc()

    def create_interface(self):
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        title_label = ttk.Label(control_frame, text="Fixed Batch Cell Analysis", font=('Arial', 14, 'bold'))
        title_label.pack(pady=(0, 20))
        
        self.create_image_selection(control_frame)
        self.create_cell_selection(control_frame)
        self.create_transform_options(control_frame)
        self.create_action_buttons(control_frame)
        
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.create_3d_plot(plot_frame)

    def create_image_selection(self, parent):
        image_frame = ttk.LabelFrame(parent, text="Image Selection", padding=10)
        image_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(image_frame, text="Select Image:").pack(anchor=tk.W)
        self.image_combo = ttk.Combobox(image_frame, textvariable=self.selected_image, values=[], state="readonly", width=20)
        self.image_combo.pack(fill=tk.X, pady=5)
        self.image_combo.bind('<<ComboboxSelected>>', self.on_image_changed)

    def create_cell_selection(self, parent):
        cell_frame = ttk.LabelFrame(parent, text="Cell Selection", padding=10)
        cell_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Label(cell_frame, text="Select Cell:").pack(anchor=tk.W)
        self.cell_combo = ttk.Combobox(cell_frame, textvariable=self.selected_cell, values=[], state="readonly", width=15)
        self.cell_combo.pack(fill=tk.X, pady=5)
        self.cell_combo.bind('<<ComboboxSelected>>', self.on_cell_changed)
        self.cell_info_label = ttk.Label(cell_frame, text="", font=('Arial', 9))
        self.cell_info_label.pack(anchor=tk.W, pady=(5, 0))

    def create_transform_options(self, parent):
        transform_frame = ttk.LabelFrame(parent, text="Transform Options", padding=10)
        transform_frame.pack(fill=tk.X, pady=(0, 15))
        
        y_flip_check = ttk.Checkbutton(transform_frame, text="Y-axis Flip", variable=self.use_y_flip, command=self.on_transform_change)
        y_flip_check.pack(anchor=tk.W)
        
        # Add cell outline option
        outline_check = ttk.Checkbutton(transform_frame, text="Show Cell Outline", variable=self.show_cell_outline, command=self.on_outline_toggle)
        outline_check.pack(anchor=tk.W)
        
        # Add skeleton display mode option
        skeleton_lines_check = ttk.Checkbutton(transform_frame, text="Skeleton as Lines", variable=self.skeleton_as_lines, command=self.on_skeleton_mode_toggle)
        skeleton_lines_check.pack(anchor=tk.W)
        
        # Add separator
        ttk.Separator(transform_frame, orient='horizontal').pack(fill=tk.X, pady=(10, 10))
        
        # Batch analysis options
        ttk.Label(transform_frame, text="Batch Analysis Options:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        auto_compare_check = ttk.Checkbutton(transform_frame, text="Auto Y-flip comparison", variable=self.auto_compare_yflip)
        auto_compare_check.pack(anchor=tk.W)
        auto_translate_check = ttk.Checkbutton(transform_frame, text="Use skeleton auto-translation", variable=self.auto_translate_skeleton)
        auto_translate_check.pack(anchor=tk.W)

    

    def create_action_buttons(self, parent):
        button_frame = ttk.LabelFrame(parent, text="Actions", padding=10)
        button_frame.pack(fill=tk.X, pady=(0, 15))
        
        calculate_btn = ttk.Button(button_frame, text="Calculate Distance", command=self.calculate_and_show_results)
        calculate_btn.pack(fill=tk.X, pady=5)
        
        batch_btn = ttk.Button(button_frame, text="Batch Distance Analysis", command=self.batch_distance_analysis)
        batch_btn.pack(fill=tk.X, pady=5)
        
        outlier_btn = ttk.Button(button_frame, text="Analyse Outliers", command=self.analyse_outliers)
        outlier_btn.pack(fill=tk.X, pady=5)

    def create_3d_plot(self, parent):
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True) 

    def on_image_changed(self, event=None):
        image_name = self.selected_image.get()
        if image_name and image_name in self.available_images:
            self.current_image = image_name
            cell_names = list(self.available_images[image_name].keys())
            self.cell_combo['values'] = cell_names
            if cell_names:
                self.selected_cell.set(cell_names[0])
                self.on_cell_changed()

    def on_cell_changed(self, event=None):
        image_name = self.selected_image.get()
        cell_name = self.selected_cell.get()
        if image_name and cell_name and image_name in self.available_images and cell_name in self.available_images[image_name]:
            self.load_cell_data_analyze_method(image_name, cell_name)
            self.update_cell_info()
            self.update_3d_plot()

    def on_transform_change(self):
        if self.current_image and self.current_cell:
            self.load_cell_data_analyze_method(self.current_image, self.current_cell)
            self.current_skeleton_translation = np.array([0, 0, 0])
            self.update_3d_plot()

    def on_outline_toggle(self):
        """Handle cell outline visibility toggle"""
        self.update_3d_plot()

    def on_skeleton_mode_toggle(self):
        """Handle skeleton display mode toggle"""
        self.update_3d_plot()

    def update_cell_info(self):
        if self.current_cell and hasattr(self, 'cell_info_label'):
            spots_count = len(self.current_spots) if self.current_spots is not None else 0
            skeleton_count = len(self.current_skeleton) if self.current_skeleton is not None else 0
            outline_count = len(self.current_outline) if self.current_outline is not None else 0
            info_text = f"Spots: {spots_count}, Skeleton: {skeleton_count}, Outline: {outline_count}"
            self.cell_info_label.config(text=info_text)


    def load_skeleton_from_vtk(self, skeleton_dir, cell_name):
        """
        Load skeleton polylines directly from VTK file using native VTK.
        This preserves the correct topology including loops and branches.
        """
        cell_base = cell_name.replace('Cell_', 'cell_').replace('cell_', '').zfill(3)
        
        # Find VTK file
        vtk_file = None
        possible_names = [
            f'cell_{cell_base}_skeleton.vtk',  # Standard format: cell_001_skeleton.vtk
            f'cell_{cell_base}.vtk',           # Alternative: cell_001.vtk
        ]
        
        for vtk_name in possible_names:
            vtk_path = os.path.join(skeleton_dir, vtk_name)
            if os.path.exists(vtk_path):
                vtk_file = vtk_path
                break
        
        # If not found, search for any VTK file matching the pattern
        if not vtk_file:
            for filename in os.listdir(skeleton_dir):
                if filename.startswith(f'cell_{cell_base}') and filename.endswith('.vtk'):
                    vtk_file = os.path.join(skeleton_dir, filename)
                    break
        
        if not vtk_file:
            print(f"Warning: VTK file not found for {cell_name}")
            return None
        
        # Read VTK file using native VTK
        polylines = self.extract_polylines_from_vtk(vtk_file)
        if polylines is None:
            return None
        
        return polylines

    def extract_polylines_from_vtk(self, vtk_file):
        """
        Extract polylines from VTK file using native VTK library.
        Returns list of polylines, each as an array of 3D points.
        """
        # Read VTK file
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(vtk_file)
        reader.Update()
        
        polydata = reader.GetOutput()
        points = polydata.GetPoints()
        lines = polydata.GetLines()
        
        if points is None or lines is None:
            print(f"Warning: No valid polyline data found in {vtk_file}")
            return None
        
        # Extract all points
        points_array = []
        for i in range(points.GetNumberOfPoints()):
            points_array.append(np.array(points.GetPoint(i)))
        
        # Extract polylines
        lines.InitTraversal()
        id_list = vtk.vtkIdList()
        
        polylines = []
        while lines.GetNextCell(id_list):
            line_points = []
            for i in range(id_list.GetNumberOfIds()):
                point_id = id_list.GetId(i)
                line_points.append(points_array[point_id])
            if len(line_points) > 1:  # Only include lines with multiple points
                polylines.append(np.array(line_points))
        
        if len(polylines) == 0:
            print(f"Warning: No polylines extracted from {vtk_file}")
            return None
        
        return polylines

    def apply_coordinate_transform_to_polylines(self, vtk_polylines, mapping_data=None, cell_name="cell_001", pixel_size_xy=0.0645):
        """
        Apply coordinate transformation to VTK polylines to match spots coordinate system.
        """
        if vtk_polylines is None:
            return None
        
        try:
            transformed_polylines = []
            
            for polyline in vtk_polylines:
                if len(polyline) == 0:
                    continue
                
                # Apply coordinate transformation (same logic as skeleton points)
                transformed_polyline = polyline.copy()
                
                # Apply mapping offset if available
                if mapping_data and cell_name in mapping_data:
                    crop_info = mapping_data[cell_name]['crop_region']
                    x_offset = crop_info['x_offset']
                    y_offset = crop_info['y_offset']
                    
                    offset_x_um = x_offset * pixel_size_xy
                    offset_y_um = y_offset * pixel_size_xy
                    
                    transformed_polyline[:, 0] += offset_x_um  # X offset
                    transformed_polyline[:, 1] += offset_y_um  # Y offset
                    # Z coordinate doesn't need offset
                
                transformed_polylines.append(transformed_polyline)
            
            return transformed_polylines
            
        except Exception as e:
            print(f"Error applying coordinate transform to polylines: {e}")
            return vtk_polylines



    def update_3d_plot(self):
        self.ax.clear()
        
        if self.current_spots is not None and len(self.current_spots) > 0:
            self.ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], self.current_spots[:, 2], 
                          c='blue', s=30, alpha=0.6, label=f'Spots ({len(self.current_spots)})')
        
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            if self.skeleton_as_lines.get() and self.skeleton_segments is not None and len(self.skeleton_segments) > 0:
                # Use VTK polylines (correct topology)
                for i, polyline in enumerate(self.skeleton_segments):
                    if len(polyline) > 1:
                        label = f'Skeleton Polylines ({len(self.skeleton_segments)})' if i == 0 else ""
                        self.ax.plot(polyline[:, 0], polyline[:, 1], polyline[:, 2], 
                                    c='red', linewidth=4, alpha=0.8, label=label)
            else:
                # Show skeleton as scatter points
                self.ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], self.current_skeleton[:, 2], 
                              c='red', s=10, alpha=0.8, label=f'Skeleton ({len(self.current_skeleton)})')
        
        if self.current_outline is not None and len(self.current_outline) > 0 and self.show_cell_outline.get():
            self.ax.plot(self.current_outline[:, 0], self.current_outline[:, 1], self.current_outline[:, 2], 
                         c='black', linewidth=2, alpha=0.7, label=f'Cell Outline ({len(self.current_outline)})')
        
        self.ax.set_xlabel('X (μm)')
        self.ax.set_ylabel('Y (μm)')
        self.ax.set_zlabel('Z (μm)')
        
        cell_title = self.current_cell if self.current_cell else "No Cell Selected"
        image_title = self.current_image if self.current_image else "No Image Selected"
        flip_status = ""
        if self.use_y_flip.get():
            flip_status += " (Y-flipped)"
        if self.use_z_flip.get():
            flip_status += " (Z-flipped)"
        
        self.ax.set_title(f'{self.channel.upper()} - {image_title} - {cell_title} - 3D Visualization{flip_status}')
        
        if (self.current_spots is not None and len(self.current_spots) > 0) or (self.current_skeleton is not None and len(self.current_skeleton) > 0) or (self.current_outline is not None and len(self.current_outline) > 0 and self.show_cell_outline.get()):
            self.ax.legend()
        
        self.set_equal_aspect_3d()
        self.canvas.draw()
        self.update_cell_info()

    def set_equal_aspect_3d(self):
        all_points = []
        
        if self.current_spots is not None and len(self.current_spots) > 0:
            all_points.append(self.current_spots)
        
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            all_points.append(self.current_skeleton)
            
        if self.current_outline is not None and len(self.current_outline) > 0 and self.show_cell_outline.get():
            all_points.append(self.current_outline)
        
        if len(all_points) > 0:
            all_points = np.vstack(all_points)
            x_range = all_points[:, 0].max() - all_points[:, 0].min()
            y_range = all_points[:, 1].max() - all_points[:, 1].min()
            z_range = all_points[:, 2].max() - all_points[:, 2].min()
            max_range = max(x_range, y_range, z_range)
            x_center = (all_points[:, 0].max() + all_points[:, 0].min()) / 2
            y_center = (all_points[:, 1].max() + all_points[:, 1].min()) / 2
            z_center = (all_points[:, 2].max() + all_points[:, 2].min()) / 2
            self.ax.set_xlim(x_center - max_range/2, x_center + max_range/2)
            self.ax.set_ylim(y_center - max_range/2, y_center + max_range/2)
            self.ax.set_zlim(z_center - max_range/2, z_center + max_range/2)


    def calculate_distances_analyze_method(self):
        """
        Calculate distances using all skeleton points (no sampling) for accurate distance calculation
        Use center-uncorrected data for visualization consistency
        """
        if not hasattr(self, 'spots_before_correction') or self.current_skeleton is None:
            return np.array([])
        
        spots_coords_um = self.spots_before_correction
        skeleton_coords = self.current_skeleton
        
        if len(spots_coords_um) == 0 or len(skeleton_coords) == 0:
            return np.array([])
        
        # Use all skeleton points for accurate distance calculation
        # Only sample spots if too many to avoid memory issues
        if len(spots_coords_um) > 500:
            spots_sample = spots_coords_um[::len(spots_coords_um)//500]
        else:
            spots_sample = spots_coords_um
        
        print(f"Distance calculation: Using {len(spots_sample)} spots and {len(skeleton_coords)} skeleton points")
        
        distances = cdist(spots_sample, skeleton_coords)
        min_distances = np.min(distances, axis=1)
        return min_distances

    def calculate_final_distances(self):
        """
        Calculate final alignment statistics using all skeleton points for accuracy
        Use center-corrected data for final statistics
        """
        if not hasattr(self, 'spots_after_correction') or self.current_skeleton is None:
            return np.array([])
        
        spots_corrected = self.spots_after_correction
        skeleton_coords = self.current_skeleton
        
        if len(spots_corrected) == 0 or len(skeleton_coords) == 0:
            return np.array([])
        
        # Use all skeleton points for accurate distance calculation
        # Only sample spots if too many to avoid memory issues
        spots_sample = spots_corrected[:500] if len(spots_corrected) > 500 else spots_corrected
        
        print(f"Final distance calculation: Using {len(spots_sample)} spots and {len(skeleton_coords)} skeleton points")
        
        distances = cdist(spots_sample, skeleton_coords)
        min_distances = np.min(distances, axis=1)
        return min_distances

    def calculate_and_show_results(self):
        if self.current_spots is None or self.current_skeleton is None:
            messagebox.showwarning("Warning", "Please select a cell first")
            return
        
        self.distances = self.calculate_distances_analyze_method()
        self.final_distances = self.calculate_final_distances()
        
        if len(self.distances) == 0:
            messagebox.showwarning("Warning", "Cannot calculate distances, please check data")
            return
            
        self.create_result_window()
        
        print(f"{self.current_image} - {self.current_cell} calculated:")
        print(f"  Graph distance (pre-correction):")
        print(f"    Mean distance: {np.mean(self.distances):.3f} μm")
        print(f"    Median distance: {np.median(self.distances):.3f} μm")
        print(f"    Standard deviation: {np.std(self.distances):.3f} μm")
        
        if len(self.final_distances) > 0:
            print(f"  Final alignment statistics (post-correction):")
            print(f"    Mean distance: {np.mean(self.final_distances):.3f} μm")
            print(f"    Median distance: {np.median(self.final_distances):.3f} μm")
            print(f"    Standard deviation: {np.std(self.final_distances):.3f} μm")

    def create_result_window(self):
        result_window = tk.Toplevel(self.root)
        result_window.title(f"Analysis Results - {self.current_image} - {self.current_cell}")
        result_window.geometry("1200x800")
        
        result_notebook = ttk.Notebook(result_window)
        result_notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        xy_frame = ttk.Frame(result_notebook)
        result_notebook.add(xy_frame, text="XY Projection")
        
        dist_frame = ttk.Frame(result_notebook)
        result_notebook.add(dist_frame, text="Distance Distribution")
        
        self.create_xy_projection_in_window(xy_frame)
        self.create_distance_distribution_in_window(dist_frame)
        
        result_notebook.select(0)
        result_window.lift()
        result_window.focus_force()

    def create_xy_projection_in_window(self, parent_frame):
        fig = Figure(figsize=(10, 8), dpi=100)
        ax = fig.add_subplot(111)
        
        if self.current_spots is not None and len(self.current_spots) > 0:
            scatter = ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], 
                               c=self.distances, cmap='viridis', s=50, alpha=0.7)
            fig.colorbar(scatter, ax=ax, label='Distance to skeleton (μm)')
        
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            if self.skeleton_as_lines.get() and self.skeleton_segments is not None and len(self.skeleton_segments) > 0:
                # Use VTK polylines for XY projection
                for i, polyline in enumerate(self.skeleton_segments):
                    if len(polyline) > 1:
                        label = 'Skeleton' if i == 0 else ""
                        ax.plot(polyline[:, 0], polyline[:, 1], 
                               c='red', linewidth=4, alpha=0.7, label=label)
            else:
                ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], 
                          c='red', s=10, alpha=0.5, label='Skeleton')
        
        if self.current_outline is not None and len(self.current_outline) > 0 and self.show_cell_outline.get():
            ax.plot(self.current_outline[:, 0], self.current_outline[:, 1], 
                    c='black', linewidth=2, alpha=0.7, label='Cell Outline')
        
        ax.set_xlabel('X (μm)')
        ax.set_ylabel('Y (μm)')
        
        cell_title = self.current_cell if self.current_cell else "Unknown Cell"
        image_title = self.current_image if self.current_image else "Unknown Image"
        flip_status = ""
        if self.use_y_flip.get():
            flip_status += " (Y-flipped)"
        if self.use_z_flip.get():
            flip_status += " (Z-flipped)"
        
        ax.set_title(f'{self.channel.upper()} - {image_title} - {cell_title} - XY Projection{flip_status}')
        
        if ((self.current_skeleton is not None and len(self.current_skeleton) > 0) or 
            (self.current_outline is not None and len(self.current_outline) > 0 and self.show_cell_outline.get())):
            ax.legend()
        
        ax.grid(True, alpha=0.3)
        ax.axis('equal')
        
        canvas = FigureCanvasTkAgg(fig, parent_frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas.draw()

    def create_distance_distribution_in_window(self, parent_frame):
        fig = Figure(figsize=(10, 8), dpi=100)
        ax = fig.add_subplot(111)
        
        if len(self.distances) > 0:
            ax.hist(self.distances, bins=20, alpha=0.7, edgecolor='black', color='skyblue')
            ax.axvline(np.mean(self.distances), color='red', linestyle='--', 
                      label=f'Mean: {np.mean(self.distances):.3f} μm')
            ax.axvline(np.median(self.distances), color='orange', linestyle='--', 
                      label=f'Median: {np.median(self.distances):.3f} μm')
        
        ax.set_xlabel('Distance to skeleton (μm)')
        ax.set_ylabel('Number of molecules')
        
        cell_title = self.current_cell if self.current_cell else "Unknown Cell"
        image_title = self.current_image if self.current_image else "Unknown Image"
        flip_status = ""
        if self.use_y_flip.get():
            flip_status += " (Y-flipped)"
        if self.use_z_flip.get():
            flip_status += " (Z-flipped)"
        
        ax.set_title(f'{self.channel.upper()} - {image_title} - {cell_title} - Distance Distribution{flip_status}')
        
        if len(self.distances) > 0:
            ax.legend()
            
            # Display chart data statistics (pre-correction, matches analyze_alignment.py chart)
            stats_text = f"""Graph Statistics (Pre-Correction):
    Count: {len(self.distances)}
    Mean: {np.mean(self.distances):.3f} μm
    Median: {np.median(self.distances):.3f} μm
    Std: {np.std(self.distances):.3f} μm
    Min: {np.min(self.distances):.3f} μm
    Max: {np.max(self.distances):.3f} μm"""

            # If there are final statistics, also display
            if hasattr(self, 'final_distances') and len(self.final_distances) > 0:
                stats_text += f"""

Final Statistics (Post-Correction):
Count: {len(self.final_distances)}
Mean: {np.mean(self.final_distances):.3f} μm
Median: {np.median(self.final_distances):.3f} μm
Std: {np.std(self.final_distances):.3f} μm"""

            stats_text += f"""

Transform: {flip_status if flip_status else 'Original'}
Note: Uses ALL skeleton points (no sampling) for accurate distance calculation"""
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)
        
        ax.grid(True, alpha=0.3)
        
        canvas = FigureCanvasTkAgg(fig, parent_frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas.draw()

    def optimize_skeleton_translation(self, spots, skeleton, search_range=2.0, step_size=0.2):
        """
        自动寻找最优骨架平移位置以最小化与spots的中位数距离
        """
        if len(spots) == 0 or len(skeleton) == 0:
            return {'optimal_skeleton': skeleton, 'translation': np.array([0, 0, 0]), 'median_distance': float('inf')}
        
        from scipy.optimize import minimize
        
        def objective_function(translation):
            """目标函数：返回平移后骨架到spots的中位数距离"""
            translated_skeleton = skeleton + translation.reshape(1, 3)
            distances = cdist(spots, translated_skeleton)
            min_distances = np.min(distances, axis=1)
            return np.median(min_distances)
        
        # 初始猜测：零平移
        initial_translation = np.array([0.0, 0.0, 0.0])
        
        # 设置搜索边界
        bounds = [(-search_range, search_range) for _ in range(3)]
        
        # 使用L-BFGS-B算法进行优化
        try:
            result = minimize(objective_function, initial_translation, method='L-BFGS-B', bounds=bounds)
            optimal_translation = result.x
            optimal_skeleton = skeleton + optimal_translation.reshape(1, 3)
            optimal_median_distance = result.fun
            
            return {
                'optimal_skeleton': optimal_skeleton,
                'translation': optimal_translation,
                'median_distance': optimal_median_distance
            }
        except Exception as e:
            # 优化失败，静默使用原始skeleton
            distances = cdist(spots, skeleton)
            min_distances = np.min(distances, axis=1)
            median_distance = np.median(min_distances)
            return {
                'optimal_skeleton': skeleton,
                'translation': np.array([0, 0, 0]),
                'median_distance': median_distance
            }


    def batch_distance_analysis(self):
        """Batch analyze distance distribution of all images and cells using the fixed correct coordinate transformation method"""
        threshold = simpledialog.askfloat("Threshold Input", "Please enter distance threshold (μm):", minvalue=0.0, initialvalue=0.5)
        if threshold is None:
            return
        
        # 创建进度条窗口
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Batch Distance Analysis")
        progress_label = ttk.Label(progress_win, text="Analyzing...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        
        all_distances = []
        summary_rows = []
        cell_spot_counts = []
        all_spots_details = []
        
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        
        for image_name, cells in self.available_images.items():
            for cell_name, cell_info in cells.items():
                try:
                    # 获取细胞编号和映射数据
                    cell_number = int(cell_name.split("_")[1])
                    cell_id_str = f"{cell_number:03d}"
                    mapping_cell_name = f'cell_{cell_id_str}'
                    mapping_data = self.coordinate_mappings.get(image_name, {})
                    
                    # 使用正确的方法加载数据
                    if self.auto_compare_yflip.get():
                        # 自动比较Y翻转模式
                        best_result = None
                        best_mean = None
                        original_y_flip = self.use_y_flip.get()
                        
                        for y_flip in [True, False]:
                            self.use_y_flip.set(y_flip)
                            
                            # 加载spots数据
                            spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                                cell_info['spots_file'], 
                                cell_number=cell_number, 
                                flip_y=y_flip, 
                                mapping_data=mapping_data,
                                silent=True
                            )
                            
                            # 加载skeleton数据
                            skeleton_coords = self.load_skeleton_txt_analyze_method(
                                cell_info['skeleton_file'], 
                                mapping_data=mapping_data,
                                cell_name=mapping_cell_name,
                                pixel_size_xy=pixel_xy/1000,
                                silent=True
                            )
                            
                            # 坐标转换
                            spots_nm_xyz = spots_coords[:, [1, 0, 2]]
                            spots_um_xyz = spots_nm_xyz / 1000.0
                            
                            # 中心校正
                            skeleton_center = np.mean(skeleton_coords, axis=0)
                            spots_center = np.mean(spots_um_xyz, axis=0)
                            center_diff = skeleton_center - spots_center
                            spots_corrected = spots_um_xyz - center_diff
                            
                            spots = spots_um_xyz  # 使用中心校正前的数据
                            skeleton = skeleton_coords
                            
                            # 应用自动平移优化（如果启用）
                            skeleton_translation = np.array([0, 0, 0])
                            if self.auto_translate_skeleton.get():
                                translation_result = self.optimize_skeleton_translation(spots, skeleton)
                                skeleton = translation_result['optimal_skeleton']
                                skeleton_translation = translation_result['translation']
                            
                            # 计算距离（使用所有skeleton点）
                            spots_sample = spots[:500] if len(spots) > 500 else spots
                            distances = cdist(spots_sample, skeleton)
                            min_distances = np.min(distances, axis=1)
                            mean_dist = np.mean(min_distances)
                            
                            if (best_mean is None) or (mean_dist < best_mean):
                                best_mean = mean_dist
                                best_result = {
                                    'spots': spots.copy(),
                                    'skeleton': skeleton.copy(),
                                    'spots_corrected': spots_corrected.copy(),
                                    'min_distances': min_distances.copy(),
                                    'y_flip': y_flip,
                                    'skeleton_translation': skeleton_translation.copy()
                                }
                        
                        # 使用最优结果
                        spots = best_result['spots']
                        skeleton = best_result['skeleton']
                        spots_corrected = best_result['spots_corrected']
                        min_distances = best_result['min_distances']
                        y_flip_used = best_result['y_flip']
                        skeleton_translation = best_result['skeleton_translation']
                        
                        # 恢复原始设置
                        self.use_y_flip.set(original_y_flip)
                    else:
                        # 使用当前界面设置
                        spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                            cell_info['spots_file'], 
                            cell_number=cell_number, 
                            flip_y=self.use_y_flip.get(), 
                            mapping_data=mapping_data,
                            silent=True
                        )
                        
                        skeleton_coords = self.load_skeleton_txt_analyze_method(
                            cell_info['skeleton_file'], 
                            mapping_data=mapping_data,
                            cell_name=mapping_cell_name,
                            pixel_size_xy=pixel_xy/1000,
                            silent=True
                        )
                        
                        spots_nm_xyz = spots_coords[:, [1, 0, 2]]
                        spots_um_xyz = spots_nm_xyz / 1000.0
                        
                        skeleton_center = np.mean(skeleton_coords, axis=0)
                        spots_center = np.mean(spots_um_xyz, axis=0)
                        center_diff = skeleton_center - spots_center
                        spots_corrected = spots_um_xyz - center_diff
                        
                        spots = spots_um_xyz
                        skeleton = skeleton_coords
                        y_flip_used = self.use_y_flip.get()
                        
                        skeleton_translation = np.array([0, 0, 0])
                        if self.auto_translate_skeleton.get():
                            translation_result = self.optimize_skeleton_translation(spots, skeleton)
                            skeleton = translation_result['optimal_skeleton']
                            skeleton_translation = translation_result['translation']
                        
                        spots_sample = spots[:500] if len(spots) > 500 else spots
                        distances = cdist(spots_sample, skeleton)
                        min_distances = np.min(distances, axis=1)
                    
                    # 统计结果
                    below = (min_distances < threshold)
                    for i, d in enumerate(min_distances):
                        all_distances.append(d)
                        all_spots_details.append({
                            'image': image_name,
                            'cell': cell_name,
                            'spot_index': i,
                            'spot_x_um': spots[i, 0] if i < len(spots) else 0,
                            'spot_y_um': spots[i, 1] if i < len(spots) else 0,
                            'spot_z_um': spots[i, 2] if i < len(spots) else 0,
                            'distance_to_skeleton_um': d,
                            'exceeds_threshold': d >= threshold,
                            'threshold_um': threshold,
                            'y_flip_used': y_flip_used,
                            'auto_compare_enabled': self.auto_compare_yflip.get(),
                            'auto_translate_enabled': self.auto_translate_skeleton.get(),
                            'skeleton_translation_x': skeleton_translation[0],
                            'skeleton_translation_y': skeleton_translation[1],
                            'skeleton_translation_z': skeleton_translation[2]
                        })
                    
                    cell_spot_counts.append(len(min_distances))
                    below_ratio = np.sum(below) / len(min_distances)
                    summary_rows.append({
                        'image': image_name,
                        'cell': cell_name,
                        'num_spots': len(min_distances),
                        'mean_distance': np.mean(min_distances),
                        'median_distance': np.median(min_distances),
                        'std_distance': np.std(min_distances),
                        'min_distance': np.min(min_distances),
                        'max_distance': np.max(min_distances),
                        f'ratio_below_{threshold}': below_ratio,
                        'y_flip_used': y_flip_used,
                        'auto_compare_enabled': self.auto_compare_yflip.get(),
                        'auto_translate_enabled': self.auto_translate_skeleton.get(),
                        'skeleton_translation_x': skeleton_translation[0],
                        'skeleton_translation_y': skeleton_translation[1],
                        'skeleton_translation_z': skeleton_translation[2]
                    })
                    
                except Exception as e:
                    # 静默跳过无法处理的细胞（例如未找到spots数据）
                    pass
                
                done += 1
                progress_var.set(done / total)
                self.root.update()
        
        # 保存结果
        output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
        os.makedirs(output_dir, exist_ok=True)
        
        dist_df = pd.DataFrame({'distance': all_distances})
        dist_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_distances_{self.skeleton_type}.csv'), index=False)
        
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_summary_{self.skeleton_type}.csv'), index=False)
        
        if all_spots_details:
            spots_details_df = pd.DataFrame(all_spots_details)
            spots_details_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_spots_details_{self.skeleton_type}.csv'), index=False)
            
            exceeds_threshold_df = spots_details_df[spots_details_df['exceeds_threshold'] == True]
            exceeds_threshold_df.to_csv(os.path.join(output_dir, f'{self.channel}_spots_exceed_threshold_{threshold}_{self.skeleton_type}.csv'), index=False)
        
        # 绘制聚合直方图和箱线图
        if not dist_df.empty:
            mean_val = dist_df['distance'].mean()
            median_val = dist_df['distance'].median()
            below_ratio = (dist_df['distance'] < threshold).mean()
            cell_count = len(cell_spot_counts)
            mean_spots_per_cell = np.mean(cell_spot_counts) if cell_spot_counts else 0
            
            plt.figure(figsize=(10, 6))
            n, bins, patches = plt.hist(dist_df['distance'], bins=30, color='skyblue', edgecolor='black', alpha=0.7)
            plt.xlabel('Distance to skeleton (μm)')
            plt.ylabel('Number of molecules')
            plt.title(f'{self.channel.upper()} - All Cells Distance Distribution (N={len(dist_df)})')
            plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
            plt.axvline(median_val, color='orange', linestyle='--', label=f'Median: {median_val:.2f}')
            plt.axvline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            plt.text(0.98, 0.95, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().transAxes, ha='right', va='top', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_histogram_{self.skeleton_type}.png'))
            plt.close()
            
            plt.figure(figsize=(8, 6))
            plt.boxplot(dist_df['distance'], vert=True, patch_artist=True, boxprops=dict(facecolor='lightgreen'))
            plt.ylabel('Distance to skeleton (μm)')
            plt.title(f'{self.channel.upper()} - All Cells Distance Boxplot ({self.skeleton_type}) (N={len(dist_df)})')
            plt.axhline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            plt.text(1.05, threshold, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().get_yaxis_transform(which='grid'), ha='left', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_boxplot_{self.skeleton_type}.png'))
            plt.close()
        
        progress_win.destroy()
        exceed_count = len(exceeds_threshold_df) if 'exceeds_threshold_df' in locals() else 0
        total_spots = len(all_spots_details) if all_spots_details else 0
        messagebox.showinfo("Batch Analysis Completed", f"Distance distribution histogram and boxplot for all cells saved to:\n{output_dir}\n\nTotal spots: {total_spots}\nSpots exceeding threshold {threshold} μm: {exceed_count}")

    def analyse_outliers(self):
        """Analyze outliers and display images and cells containing outliers in a new window"""
        # Ask for threshold
        threshold = simpledialog.askfloat("Outlier Threshold Input", "Please enter outlier distance threshold (μm):", minvalue=0.0, initialvalue=0.5)
        if threshold is None:
            return
        
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Analyse Outliers Progress")
        progress_label = ttk.Label(progress_win, text="Analyzing outliers...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        
        outlier_images = {}  # {image_name: {cell_name: outlier_info}}
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        
        for image_name, cells in self.available_images.items():
            for cell_name, cell_info in cells.items():
                try:
                    cell_number = int(cell_name.split("_")[1])
                    cell_id_str = f"{cell_number:03d}"
                    mapping_cell_name = f'cell_{cell_id_str}'
                    mapping_data = self.coordinate_mappings.get(image_name, {})
                    
                    spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                        cell_info['spots_file'], 
                        cell_number=cell_number, 
                        flip_y=self.use_y_flip.get(), 
                        mapping_data=mapping_data,
                        silent=True
                    )
                    
                    skeleton_coords = self.load_skeleton_txt_analyze_method(
                        cell_info['skeleton_file'], 
                        mapping_data=mapping_data,
                        cell_name=mapping_cell_name,
                        pixel_size_xy=pixel_xy/1000,
                        silent=True
                    )
                    
                    spots_nm_xyz = spots_coords[:, [1, 0, 2]]
                    spots_um_xyz = spots_nm_xyz / 1000.0
                    
                    skeleton_center = np.mean(skeleton_coords, axis=0)
                    spots_center = np.mean(spots_um_xyz, axis=0)
                    center_diff = skeleton_center - spots_center
                    spots_corrected = spots_um_xyz - center_diff
                    
                    spots = spots_um_xyz 
                    skeleton = skeleton_coords
                    
                    skeleton_translation = np.array([0, 0, 0])
                    if self.auto_translate_skeleton.get():
                        translation_result = self.optimize_skeleton_translation(spots, skeleton)
                        skeleton = translation_result['optimal_skeleton']
                        skeleton_translation = translation_result['translation']
                    
                    # Load VTK polylines for skeleton
                    skeleton_dir = os.path.dirname(cell_info['skeleton_file'])
                    vtk_polylines = self.load_skeleton_from_vtk(skeleton_dir, cell_name)
                    skeleton_segments = self.apply_coordinate_transform_to_polylines(
                        vtk_polylines,
                        mapping_data=mapping_data,
                        cell_name=mapping_cell_name,
                        pixel_size_xy=pixel_xy/1000
                    )
                    
                    distances = cdist(spots, skeleton)
                    min_distances = np.min(distances, axis=1)
                    
                    outlier_indices = np.where(min_distances >= threshold)[0]
                    if len(outlier_indices) > 0:
                        outlier_images[image_name] = outlier_images.get(image_name, {})
                        outlier_images[image_name][cell_name] = {
                            'spots': spots,
                            'skeleton': skeleton,
                            'skeleton_segments': skeleton_segments,
                            'distances': min_distances,
                            'outlier_indices': outlier_indices,
                            'threshold': threshold,
                            'skeleton_translation': skeleton_translation,
                            'original_cell_info': cell_info
                        }
                
                except Exception as e:
                    pass
                
                done += 1
                progress_var.set(done / total)
                self.root.update()
        
        progress_win.destroy()
        
        if not outlier_images:
            messagebox.showinfo("Analysis Results", f"No outliers found exceeding threshold {threshold} μm")
            return
        
        total_outlier_cells = sum(len(cells) for cells in outlier_images.values())
        total_outliers = sum(len(cell_data['outlier_indices']) for image_cells in outlier_images.values() 
                           for cell_data in image_cells.values())
        
        self.create_outlier_analysis_window(outlier_images, threshold, total_outlier_cells, total_outliers)

    def create_outlier_analysis_window(self, outlier_images, threshold, total_outlier_cells, total_outliers):
        outlier_window = tk.Toplevel(self.root)
        outlier_window.title(f"Outlier Analysis (Fixed Method) - Threshold: {threshold} μm")
        outlier_window.geometry("1500x950")
        
        main_frame = ttk.Frame(outlier_window)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        title_label = ttk.Label(control_frame, text=f"Outlier Analysis (Fixed)\nThreshold: {threshold} μm", 
                               font=('Arial', 14, 'bold'))
        title_label.pack(pady=(0, 20))
        
        stats_frame = ttk.LabelFrame(control_frame, text="Statistics", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 15))
        
        stats_text = f"Images with outliers: {len(outlier_images)}\nCells with outliers: {total_outlier_cells}\nTotal outlier spots: {total_outliers}"
        stats_label = ttk.Label(stats_frame, text=stats_text, font=('Arial', 10))
        stats_label.pack(anchor=tk.W)
        
        image_frame = ttk.LabelFrame(control_frame, text="Image Selection", padding=10)
        image_frame.pack(fill=tk.X, pady=(0, 15))
        
        outlier_selected_image = tk.StringVar()
        ttk.Label(image_frame, text="Select Image:").pack(anchor=tk.W)
        image_combo = ttk.Combobox(image_frame, textvariable=outlier_selected_image, 
                                  values=sorted(list(outlier_images.keys())), state="readonly", width=20)
        image_combo.pack(fill=tk.X, pady=5)
        
        cell_frame = ttk.LabelFrame(control_frame, text="Cell Selection", padding=10)
        cell_frame.pack(fill=tk.X, pady=(0, 15))
        
        outlier_selected_cell = tk.StringVar()
        ttk.Label(cell_frame, text="Select Cell:").pack(anchor=tk.W)
        cell_combo = ttk.Combobox(cell_frame, textvariable=outlier_selected_cell, 
                                 values=[], state="readonly", width=15)
        cell_combo.pack(fill=tk.X, pady=5)
        
        cell_info_label = ttk.Label(cell_frame, text="", font=('Arial', 9))
        cell_info_label.pack(anchor=tk.W, pady=(5, 0))
        
        # Add skeleton display options for outlier analysis
        skeleton_frame = ttk.LabelFrame(control_frame, text="Skeleton Display", padding=10)
        skeleton_frame.pack(fill=tk.X, pady=(0, 15))
        
        outlier_skeleton_as_lines = tk.BooleanVar(value=True)
        skeleton_lines_check = ttk.Checkbutton(skeleton_frame, text="Skeleton as Lines", variable=outlier_skeleton_as_lines)
        skeleton_lines_check.pack(anchor=tk.W)
        
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        outlier_fig = Figure(figsize=(10, 8), dpi=100)
        outlier_ax = outlier_fig.add_subplot(111, projection='3d')
        outlier_canvas = FigureCanvasTkAgg(outlier_fig, plot_frame)
        outlier_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        def update_outlier_cell_selection(event=None):
            image_name = outlier_selected_image.get()
            if image_name and image_name in outlier_images:
                cell_names = list(outlier_images[image_name].keys())
                cell_combo['values'] = cell_names
                if cell_names:
                    outlier_selected_cell.set(cell_names[0])
                    update_outlier_plot()
        
        def update_outlier_plot(event=None):
            image_name = outlier_selected_image.get()
            cell_name = outlier_selected_cell.get()
            
            if not image_name or not cell_name or image_name not in outlier_images or cell_name not in outlier_images[image_name]:
                return
            
            cell_data = outlier_images[image_name][cell_name]
            spots = cell_data['spots']
            skeleton = cell_data['skeleton']
            distances = cell_data['distances']
            outlier_indices = cell_data['outlier_indices']
            threshold_val = cell_data['threshold']
            
            outlier_ax.clear()
            
            normal_indices = np.where(distances < threshold_val)[0]
            if len(normal_indices) > 0:
                normal_spots = spots[normal_indices]
                outlier_ax.scatter(normal_spots[:, 0], normal_spots[:, 1], normal_spots[:, 2], 
                                 c='blue', s=30, alpha=0.6, label=f'Normal spots ({len(normal_indices)})')
            
            if len(outlier_indices) > 0:
                outlier_spots = spots[outlier_indices]
                outlier_ax.scatter(outlier_spots[:, 0], outlier_spots[:, 1], outlier_spots[:, 2], 
                                 c='red', s=50, alpha=0.8, label=f'Outlier spots ({len(outlier_indices)})')
            
            if len(skeleton) > 0:
                skeleton_segments = cell_data.get('skeleton_segments', None)
                if outlier_skeleton_as_lines.get() and skeleton_segments is not None and len(skeleton_segments) > 0:
                    # Use VTK polylines for skeleton display
                    for i, polyline in enumerate(skeleton_segments):
                        if len(polyline) > 1:
                            label = f'Skeleton Polylines ({len(skeleton_segments)})' if i == 0 else ""
                            outlier_ax.plot(polyline[:, 0], polyline[:, 1], polyline[:, 2], 
                                           c='gray', linewidth=4, alpha=0.7, label=label)
                else:
                    # Show skeleton as scatter points
                    outlier_ax.scatter(skeleton[:, 0], skeleton[:, 1], skeleton[:, 2], 
                                     c='gray', s=10, alpha=0.5, label=f'Skeleton ({len(skeleton)})')
            
            outlier_ax.set_xlabel('X (μm)')
            outlier_ax.set_ylabel('Y (μm)')
            outlier_ax.set_zlabel('Z (μm)')
            
            flip_status = ""
            if self.use_y_flip.get():
                flip_status += " (Y-flipped)"
            if self.use_z_flip.get():
                flip_status += " (Z-flipped)"
            
            outlier_ax.set_title(f'{self.channel.upper()} - {image_name} - {cell_name} - Outlier Analysis (Fixed){flip_status}\nThreshold: {threshold_val} μm')
            outlier_ax.legend()
            
            if len(spots) > 0:
                all_points = np.vstack([spots, skeleton])
                x_range = all_points[:, 0].max() - all_points[:, 0].min()
                y_range = all_points[:, 1].max() - all_points[:, 1].min()
                z_range = all_points[:, 2].max() - all_points[:, 2].min()
                max_range = max(x_range, y_range, z_range)
                x_center = (all_points[:, 0].max() + all_points[:, 0].min()) / 2
                y_center = (all_points[:, 1].max() + all_points[:, 1].min()) / 2
                z_center = (all_points[:, 2].max() + all_points[:, 2].min()) / 2
                outlier_ax.set_xlim(x_center - max_range/2, x_center + max_range/2)
                outlier_ax.set_ylim(y_center - max_range/2, y_center + max_range/2)
                outlier_ax.set_zlim(z_center - max_range/2, z_center + max_range/2)
            
            outlier_canvas.draw()
            
            total_spots = len(spots)
            outlier_count = len(outlier_indices)
            normal_count = len(normal_indices)
            outlier_ratio = outlier_count / total_spots * 100 if total_spots > 0 else 0
            
            info_text = f"Total: {total_spots}, Normal: {normal_count}, Outliers: {outlier_count} ({outlier_ratio:.1f}%)"
            cell_info_label.config(text=info_text)
        
        image_combo.bind('<<ComboboxSelected>>', update_outlier_cell_selection)
        cell_combo.bind('<<ComboboxSelected>>', update_outlier_plot)
        skeleton_lines_check.config(command=update_outlier_plot)

        if outlier_images:
            first_image = sorted(list(outlier_images.keys()))[0]
            outlier_selected_image.set(first_image)
            update_outlier_cell_selection()
        
        action_frame = ttk.LabelFrame(control_frame, text="Actions", padding=10)
        action_frame.pack(fill=tk.X, pady=(15, 0))
        
        def save_outlier_details():
            outlier_details = []
            for image_name, image_cells in outlier_images.items():
                for cell_name, cell_data in image_cells.items():
                    spots = cell_data['spots']
                    distances = cell_data['distances']
                    outlier_indices = cell_data['outlier_indices']
                    
                    for idx in outlier_indices:
                        outlier_details.append({
                            'image': image_name,
                            'cell': cell_name,
                            'spot_index': idx,
                            'spot_x_um': spots[idx, 0],
                            'spot_y_um': spots[idx, 1],
                            'spot_z_um': spots[idx, 2],
                            'distance_to_skeleton_um': distances[idx],
                            'threshold_um': threshold,
                            'y_flip_used': self.use_y_flip.get(),
                            'z_flip_used': self.use_z_flip.get(),
                            'auto_translate_enabled': self.auto_translate_skeleton.get(),
                            'skeleton_translation_x': cell_data['skeleton_translation'][0],
                            'skeleton_translation_y': cell_data['skeleton_translation'][1],
                            'skeleton_translation_z': cell_data['skeleton_translation'][2]
                        })
            
            if outlier_details:
                output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
                os.makedirs(output_dir, exist_ok=True)
                outlier_df = pd.DataFrame(outlier_details)
                outlier_file = os.path.join(output_dir, f'{self.channel}_outlier_analysis_{threshold}_{self.skeleton_type}.csv')
                outlier_df.to_csv(outlier_file, index=False)
                messagebox.showinfo("Save Successful", f"Outlier detailed information saved to:\n{outlier_file}")
        
        save_all_btn = ttk.Button(action_frame, text="Save All Outlier Details", command=save_outlier_details)
        save_all_btn.pack(fill=tk.X, pady=2)


def main():
    root = tk.Tk()
    app = Interactive3DBatchVisualizerFixed(root)
    root.mainloop()

if __name__ == "__main__":
    main() 