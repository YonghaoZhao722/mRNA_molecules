# Dependencies
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import json
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
import cv2
from scipy.spatial import ConvexHull

# Try to import optional dependencies
try:
    import alphashape
    ALPHASHAPE_AVAILABLE = True
except ImportError:
    print("Warning: alphashape library not available. Using convex hull as fallback.")
    ALPHASHAPE_AVAILABLE = False

class AlphaShapeProjectionAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Alpha Shape Projection Analysis")
        self.root.geometry("1600x1000")

        # Path settings (same as original file)
        self.base_dir = r'Y333 ATP6 ATP2'
        self.skeleton_type = 'extracted_cells_conn_nonadaptive_rm'
        self.skeleton_root = os.path.join(self.base_dir, self.skeleton_type)
        self.channel = 'atp6_corrected'
        self.spots_root = os.path.join(self.base_dir, f'{self.channel}_spots')
        self.mask_root = os.path.join(self.base_dir, 'aligned_masks')

        # Data structures
        self.available_images = {}
        self.coordinate_mappings = {}
        self.selected_image = tk.StringVar()
        self.selected_cell = tk.StringVar()
        self.current_cell = None
        self.current_image = None

        # Data variables
        self.current_spots = None
        self.current_skeleton = None
        self.use_y_flip = tk.BooleanVar(value=False)
        
        # Alpha shape parameters
        self.alpha_value = tk.DoubleVar(value=0.1)
        self.resolution = tk.IntVar(value=100)  # Grid resolution for mask
        self.voxel_resolution = tk.IntVar(value=50)  # 3D voxel resolution for volume
        self.use_spots_for_shape = tk.BooleanVar(value=True)
        self.use_skeleton_for_shape = tk.BooleanVar(value=True)
        
        # Results storage
        self.alpha_shapes = {}  # {'xy': shape, 'xz': shape, 'yz': shape}
        self.masks = {}  # {'xy': mask, 'xz': mask, 'yz': mask}
        self.volume_3d = None  # 3D volume reconstructed from projections
        self.volume_value = 0.0  # Computed volume in cubic micrometers
        
        # Create interface
        self.create_interface()
        # Scan available data
        self.scan_available_images_and_cells()
        # Auto load first data
        if self.available_images:
            first_image = list(self.available_images.keys())[0]
            self.selected_image.set(first_image)
            self.on_image_changed()

    def load_spots_fishquant_analyze_method(self, file_path, cell_number=1, flip_y=True, mapping_data=None, silent=False):
        """
        Copy of the correct load_spots_fishquant method from the original file
        """
        # Read file to find SPOTS section
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find pixel size information
        pixel_xy, pixel_z = None, None
        for i, line in enumerate(lines):
            if line.startswith('Pix-XY'):
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
            
            if line_stripped == target_cell:
                cell_found = True
                if not silent:
                    print(f"Found target cell: Cell_{cell_number}")
                continue
            
            if cell_found and line.startswith('Pos_Y'):
                spots_start = i
                continue
            
            if cell_found and spots_start != -1 and line.startswith('CELL'):
                spots_end = i
                break
                
            if cell_found and spots_start == -1 and line.startswith('CELL'):
                if not silent:
                    print(f"Cell_{cell_number} has no spots data")
                raise ValueError(f"Cell_{cell_number} has no spots data")
        
        if cell_found and spots_start != -1 and spots_end == -1:
            spots_end = len(lines)
        
        if not cell_found:
            raise ValueError(f"Cell_{cell_number} data not found")
        
        if spots_start == -1:
            raise ValueError(f"Cell_{cell_number} has no spots data")
        
        # Read spots data
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
        
        # Y-axis flip if needed
        if flip_y and mapping_data:
            cell_name = f"cell_{cell_number:03d}"
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                y_start = crop_info['y_start']
                y_end = crop_info['y_end']
                
                y_start_nm = y_start * pixel_xy
                y_end_nm = y_end * pixel_xy
                
                flip_center = y_start_nm + y_end_nm
                coords[:, 0] = flip_center - coords[:, 0]
                
                if not silent:
                    print(f"Y-axis flip applied")
        
        if not silent:
            print(f"Cell_{cell_number} spots loaded: {len(coords)} points")
        
        return coords, pixel_xy, pixel_z

    def load_skeleton_txt_analyze_method(self, file_path, mapping_data=None, cell_name="cell_001", pixel_size_xy=0.0645, silent=False):
        """
        Copy of the correct load_skeleton_txt method from the original file
        """
        df = pd.read_csv(file_path, sep='\t')
        coords = df[['x', 'y', 'z']].values
        
        if not silent:
            print(f"Skeleton data loaded: {len(coords)} points")
        
        # Convert to absolute coordinates if mapping data available
        if mapping_data and cell_name:
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                x_offset = crop_info['x_offset']
                y_offset = crop_info['y_offset']
                
                offset_x_um = x_offset * pixel_size_xy
                offset_y_um = y_offset * pixel_size_xy
                
                coords[:, 0] += offset_x_um
                coords[:, 1] += offset_y_um
                
                if not silent:
                    print(f"Applied coordinate offset: X+{offset_x_um:.3f}μm, Y+{offset_y_um:.3f}μm")
        
        return coords

    def scan_available_images_and_cells(self):
        """Scan available images and cells (same logic as original)"""
        self.available_images = {}
        if not os.path.exists(self.skeleton_root):
            print(f"Skeleton root not found: {self.skeleton_root}")
            return
            
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
                continue

            # Scan skeleton files
            skeleton_files = [f for f in os.listdir(image_path) if f.startswith('cell_') and f.endswith('.txt')]
            image_cells = {}
            
            for skel_file in skeleton_files:
                try:
                    cell_num = int(skel_file.split('_')[1].split('.')[0])
                    cell_name = f"Cell_{cell_num}"
                    
                    # Test if cell has spots data
                    try:
                        test_coords, _, _ = self.load_spots_fishquant_analyze_method(
                            matched_spot, 
                            cell_number=cell_num, 
                            flip_y=False,
                            mapping_data=None,
                            silent=True
                        )
                        
                        if len(test_coords) > 0:
                            image_cells[cell_name] = {
                                'spots_file': matched_spot,
                                'skeleton_file': os.path.join(image_path, skel_file)
                            }
                            print(f"Found {image_folder}-{cell_name} ({len(test_coords)} spots)")
                            
                    except ValueError:
                        pass
                        
                except Exception as e:
                    pass
            
            if image_cells:
                self.available_images[image_folder] = image_cells
        
        # Update dropdown
        if hasattr(self, 'image_combo'):
            self.image_combo['values'] = list(self.available_images.keys())
            if self.available_images:
                self.selected_image.set(list(self.available_images.keys())[0])

    def load_cell_data(self, image_name, cell_name):
        """Load cell data using the correct method from original file"""
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
            
            # Load spots data
            spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                cell_info['spots_file'], 
                cell_number=cell_number, 
                flip_y=self.use_y_flip.get(), 
                mapping_data=mapping_data
            )
            
            # Load skeleton data
            skeleton_coords = self.load_skeleton_txt_analyze_method(
                cell_info['skeleton_file'], 
                mapping_data=mapping_data,
                cell_name=mapping_cell_name,
                pixel_size_xy=pixel_xy/1000  # Convert to micrometers
            )
            
            # Coordinate transformation (same as original)
            spots_nm_xyz = spots_coords[:, [1, 0, 2]]  # from (y, x, z) to (x, y, z)
            spots_um_xyz = spots_nm_xyz / 1000.0  # Convert to micrometers
            
            # Center correction
            skeleton_center = np.mean(skeleton_coords, axis=0)
            spots_center = np.mean(spots_um_xyz, axis=0)
            center_diff = skeleton_center - spots_center
            
            # Store data
            self.current_spots = spots_um_xyz - center_diff  # Use center-corrected data
            self.current_skeleton = skeleton_coords.copy()
            
            print(f"Data loaded successfully:")
            print(f"  Spots: {len(self.current_spots)} points")
            print(f"  Skeleton: {len(self.current_skeleton)} points")
            
            # Update window title
            self.root.title(f"Alpha Shape Analysis - {image_name} - {cell_name}")
            
        except Exception as e:
            print(f"Failed to load data: {e}")
            import traceback
            traceback.print_exc()

    def create_interface(self):
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Control panel
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        title_label = ttk.Label(control_frame, text="Alpha Shape Projection Analysis", 
                               font=('Arial', 14, 'bold'))
        title_label.pack(pady=(0, 20))
        
        self.create_data_selection(control_frame)
        self.create_alpha_parameters(control_frame)
        self.create_action_buttons(control_frame)
        
        # Plot area
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.create_plot_area(plot_frame)

    def create_data_selection(self, parent):
        data_frame = ttk.LabelFrame(parent, text="Data Selection", padding=10)
        data_frame.pack(fill=tk.X, pady=(0, 15))
        
        # Image selection
        ttk.Label(data_frame, text="Select Image:").pack(anchor=tk.W)
        self.image_combo = ttk.Combobox(data_frame, textvariable=self.selected_image, 
                                       values=[], state="readonly", width=25)
        self.image_combo.pack(fill=tk.X, pady=5)
        self.image_combo.bind('<<ComboboxSelected>>', self.on_image_changed)
        
        # Cell selection
        ttk.Label(data_frame, text="Select Cell:").pack(anchor=tk.W, pady=(10,0))
        self.cell_combo = ttk.Combobox(data_frame, textvariable=self.selected_cell, 
                                      values=[], state="readonly", width=15)
        self.cell_combo.pack(fill=tk.X, pady=5)
        self.cell_combo.bind('<<ComboboxSelected>>', self.on_cell_changed)
        
        # Options
        y_flip_check = ttk.Checkbutton(data_frame, text="Y-axis Flip", 
                                      variable=self.use_y_flip, command=self.on_data_change)
        y_flip_check.pack(anchor=tk.W, pady=(10,0))

    def create_alpha_parameters(self, parent):
        alpha_frame = ttk.LabelFrame(parent, text="Alpha Shape Parameters", padding=10)
        alpha_frame.pack(fill=tk.X, pady=(0, 15))
        
        # Alpha value
        alpha_label_frame = ttk.Frame(alpha_frame)
        alpha_label_frame.pack(fill=tk.X, pady=2)
        ttk.Label(alpha_label_frame, text="Alpha Value:", width=12).pack(side=tk.LEFT)
        self.alpha_label = ttk.Label(alpha_label_frame, text="0.10", width=8)
        self.alpha_label.pack(side=tk.RIGHT)
        
        alpha_scale = ttk.Scale(alpha_frame, from_=0.01, to=2.0, variable=self.alpha_value, 
                               orient=tk.HORIZONTAL, command=self.on_alpha_change)
        alpha_scale.pack(fill=tk.X, pady=2)
        
        # Resolution
        res_label_frame = ttk.Frame(alpha_frame)
        res_label_frame.pack(fill=tk.X, pady=2)
        ttk.Label(res_label_frame, text="Resolution:", width=12).pack(side=tk.LEFT)
        self.res_label = ttk.Label(res_label_frame, text="100", width=8)
        self.res_label.pack(side=tk.RIGHT)
        
        res_scale = ttk.Scale(alpha_frame, from_=50, to=500, variable=self.resolution, 
                             orient=tk.HORIZONTAL, command=self.on_resolution_change)
        res_scale.pack(fill=tk.X, pady=2)
        
        # Voxel resolution for 3D volume
        voxel_label_frame = ttk.Frame(alpha_frame)
        voxel_label_frame.pack(fill=tk.X, pady=2)
        ttk.Label(voxel_label_frame, text="Voxel Res:", width=12).pack(side=tk.LEFT)
        self.voxel_label = ttk.Label(voxel_label_frame, text="50", width=8)
        self.voxel_label.pack(side=tk.RIGHT)
        
        voxel_scale = ttk.Scale(alpha_frame, from_=20, to=100, variable=self.voxel_resolution, 
                               orient=tk.HORIZONTAL, command=self.on_voxel_resolution_change)
        voxel_scale.pack(fill=tk.X, pady=2)
        

        
        # Performance tip
        perf_tip = ttk.Label(alpha_frame, text="Tip: Using smooth polyhedron display for optimal 3D visualization", 
                            font=('Arial', 8), foreground='gray')
        perf_tip.pack(anchor=tk.W, pady=(2, 0))
        
        # Data source selection
        ttk.Label(alpha_frame, text="Use for Alpha Shape:", font=('Arial', 9, 'bold')).pack(anchor=tk.W, pady=(10,5))
        spots_check = ttk.Checkbutton(alpha_frame, text="Use Spots", 
                                     variable=self.use_spots_for_shape, command=self.on_source_change)
        spots_check.pack(anchor=tk.W)
        skeleton_check = ttk.Checkbutton(alpha_frame, text="Use Skeleton", 
                                        variable=self.use_skeleton_for_shape, command=self.on_source_change)
        skeleton_check.pack(anchor=tk.W)

    def create_action_buttons(self, parent):
        button_frame = ttk.LabelFrame(parent, text="Actions", padding=10)
        button_frame.pack(fill=tk.X, pady=(0, 15))
        
        compute_btn = ttk.Button(button_frame, text="Compute Alpha Shapes", 
                                command=self.compute_alpha_shapes)
        compute_btn.pack(fill=tk.X, pady=2)
        
        save_shapes_btn = ttk.Button(button_frame, text="Save Alpha Shapes", 
                                    command=self.save_alpha_shapes)
        save_shapes_btn.pack(fill=tk.X, pady=2)
        
        save_masks_btn = ttk.Button(button_frame, text="Save Masks", 
                                   command=self.save_masks)
        save_masks_btn.pack(fill=tk.X, pady=2)
        
        batch_btn = ttk.Button(button_frame, text="Batch Process All Cells", 
                              command=self.batch_process)
        batch_btn.pack(fill=tk.X, pady=5)
        
        # 3D Volume analysis buttons
        ttk.Separator(button_frame, orient='horizontal').pack(fill=tk.X, pady=5)
        
        volume_btn = ttk.Button(button_frame, text="Compute 3D Polyhedron", 
                               command=self.compute_3d_volume)
        volume_btn.pack(fill=tk.X, pady=2)
        
        save_volume_btn = ttk.Button(button_frame, text="Save 3D Volume", 
                                    command=self.save_3d_volume)
        save_volume_btn.pack(fill=tk.X, pady=2)
        
        batch_volume_btn = ttk.Button(button_frame, text="Batch Volume Analysis", 
                                     command=self.batch_volume_analysis)
        batch_volume_btn.pack(fill=tk.X, pady=5)

    def create_plot_area(self, parent):
        # Create notebook for different projections
        self.notebook = ttk.Notebook(parent)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # XY projection
        self.xy_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.xy_frame, text="XY Projection")
        self.xy_fig = Figure(figsize=(8, 6), dpi=100)
        self.xy_ax = self.xy_fig.add_subplot(111)
        self.xy_canvas = FigureCanvasTkAgg(self.xy_fig, self.xy_frame)
        self.xy_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # XZ projection
        self.xz_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.xz_frame, text="XZ Projection")
        self.xz_fig = Figure(figsize=(8, 6), dpi=100)
        self.xz_ax = self.xz_fig.add_subplot(111)
        self.xz_canvas = FigureCanvasTkAgg(self.xz_fig, self.xz_frame)
        self.xz_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # YZ projection
        self.yz_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.yz_frame, text="YZ Projection")
        self.yz_fig = Figure(figsize=(8, 6), dpi=100)
        self.yz_ax = self.yz_fig.add_subplot(111)
        self.yz_canvas = FigureCanvasTkAgg(self.yz_fig, self.yz_frame)
        self.yz_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # 3D Polyhedron
        self.volume_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.volume_frame, text="3D Polyhedron")
        self.volume_fig = Figure(figsize=(8, 6), dpi=100)
        self.volume_ax = self.volume_fig.add_subplot(111, projection='3d')
        self.volume_canvas = FigureCanvasTkAgg(self.volume_fig, self.volume_frame)
        self.volume_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

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
        if (image_name and cell_name and 
            image_name in self.available_images and 
            cell_name in self.available_images[image_name]):
            self.load_cell_data(image_name, cell_name)
            self.update_plots()

    def on_data_change(self):
        if self.current_image and self.current_cell:
            self.load_cell_data(self.current_image, self.current_cell)
            self.update_plots()

    def on_alpha_change(self, value=None):
        self.alpha_label.config(text=f"{self.alpha_value.get():.2f}")
        if hasattr(self, 'alpha_shapes') and self.alpha_shapes:
            self.compute_alpha_shapes()

    def on_resolution_change(self, value=None):
        self.res_label.config(text=f"{int(self.resolution.get())}")

    def on_voxel_resolution_change(self, value=None):
        self.voxel_label.config(text=f"{int(self.voxel_resolution.get())}")

    def on_source_change(self):
        if hasattr(self, 'alpha_shapes'):
            self.compute_alpha_shapes()



    def _point_in_convex_hull(self, point, hull):
        """Check if a point is inside a convex hull"""
        try:
            # Add the point to the hull and check if the volume increases
            extended_points = np.vstack([hull.points, point.reshape(1, -1)])
            extended_hull = ConvexHull(extended_points)
            return len(extended_hull.vertices) == len(hull.vertices)
        except:
            return False

    def get_shape_data(self):
        """Get data points for alpha shape computation"""
        data_points = []
        
        if self.use_spots_for_shape.get() and self.current_spots is not None:
            data_points.append(self.current_spots)
            
        if self.use_skeleton_for_shape.get() and self.current_skeleton is not None:
            data_points.append(self.current_skeleton)
        
        if not data_points:
            return None
            
        return np.vstack(data_points)

    def compute_alpha_shapes(self):
        """Compute alpha shapes for all three projections"""
        if self.current_spots is None and self.current_skeleton is None:
            messagebox.showwarning("Warning", "No data loaded")
            return
        
        shape_data = self.get_shape_data()
        if shape_data is None:
            messagebox.showwarning("Warning", "No data source selected")
            return
        
        alpha = self.alpha_value.get()
        self.alpha_shapes = {}
        
        try:
            # XY projection (z=0)
            xy_points = shape_data[:, [0, 1]]
            if len(xy_points) >= 3:
                if ALPHASHAPE_AVAILABLE:
                    xy_alpha_shape = alphashape.alphashape(xy_points, alpha)
                else:
                    # Fallback to convex hull if alphashape is not available
                    hull_2d = ConvexHull(xy_points)
                    xy_alpha_shape = hull_2d.points[hull_2d.simplices]
                self.alpha_shapes['xy'] = xy_alpha_shape
            
            # XZ projection (y=1)  
            xz_points = shape_data[:, [0, 2]]
            if len(xz_points) >= 3:
                if ALPHASHAPE_AVAILABLE:
                    xz_alpha_shape = alphashape.alphashape(xz_points, alpha)
                else:
                    # Fallback to convex hull if alphashape is not available
                    hull_2d = ConvexHull(xz_points)
                    xz_alpha_shape = hull_2d.points[hull_2d.simplices]
                self.alpha_shapes['xz'] = xz_alpha_shape
            
            # YZ projection (x=2)
            yz_points = shape_data[:, [1, 2]]
            if len(yz_points) >= 3:
                if ALPHASHAPE_AVAILABLE:
                    yz_alpha_shape = alphashape.alphashape(yz_points, alpha)
                else:
                    # Fallback to convex hull if alphashape is not available
                    hull_2d = ConvexHull(yz_points)
                    yz_alpha_shape = hull_2d.points[hull_2d.simplices]
                self.alpha_shapes['yz'] = yz_alpha_shape
            
            print(f"Alpha shapes computed with alpha={alpha}")
            self.compute_masks()
            self.update_plots()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to compute alpha shapes: {str(e)}")
            print(f"Alpha shape computation error: {e}")

    def compute_masks(self):
        """Convert alpha shapes to binary masks"""
        if not self.alpha_shapes:
            return
        
        self.masks = {}
        resolution = int(self.resolution.get())
        
        shape_data = self.get_shape_data()
        if shape_data is None:
            return
        
        try:
            # XY mask
            if 'xy' in self.alpha_shapes and self.alpha_shapes['xy'] is not None:
                xy_points = shape_data[:, [0, 1]]
                x_min, x_max = xy_points[:, 0].min(), xy_points[:, 0].max()
                y_min, y_max = xy_points[:, 1].min(), xy_points[:, 1].max()
                
                # Add padding
                padding = 0.1
                x_range = x_max - x_min
                y_range = y_max - y_min
                x_min -= x_range * padding
                x_max += x_range * padding
                y_min -= y_range * padding
                y_max += y_range * padding
                
                x_grid = np.linspace(x_min, x_max, resolution)
                y_grid = np.linspace(y_min, y_max, resolution)
                xx, yy = np.meshgrid(x_grid, y_grid)
                
                # Create mask
                mask = np.zeros((resolution, resolution), dtype=bool)
                
                if ALPHASHAPE_AVAILABLE and hasattr(self.alpha_shapes['xy'], 'contains'):
                    # Use alphashape method
                    for i in range(resolution):
                        for j in range(resolution):
                            point = Point(xx[i, j], yy[i, j])
                            if hasattr(self.alpha_shapes['xy'], 'contains'):
                                mask[i, j] = self.alpha_shapes['xy'].contains(point)
                            elif hasattr(self.alpha_shapes['xy'], 'intersects'):
                                mask[i, j] = self.alpha_shapes['xy'].intersects(point)
                else:
                    # Fallback convex hull approach
                    hull_points = self.alpha_shapes['xy']
                    if len(hull_points) > 0:
                        # Create polygon from convex hull points
                        try:
                            polygon = Polygon(hull_points.reshape(-1, 2))
                            for i in range(resolution):
                                for j in range(resolution):
                                    point = Point(xx[i, j], yy[i, j])
                                    mask[i, j] = polygon.contains(point)
                        except:
                            # If polygon creation fails, use simple point-in-hull test
                            hull = ConvexHull(xy_points)
                            for i in range(resolution):
                                for j in range(resolution):
                                    # Simple point-in-convex-hull test
                                    test_point = np.array([xx[i, j], yy[i, j]])
                                    mask[i, j] = self._point_in_convex_hull(test_point, hull)
                
                self.masks['xy'] = {
                    'mask': mask,
                    'extent': [x_min, x_max, y_min, y_max],
                    'resolution': resolution
                }
            
            # XZ mask
            if 'xz' in self.alpha_shapes and self.alpha_shapes['xz'] is not None:
                xz_points = shape_data[:, [0, 2]]
                x_min, x_max = xz_points[:, 0].min(), xz_points[:, 0].max()
                z_min, z_max = xz_points[:, 1].min(), xz_points[:, 1].max()
                
                padding = 0.1
                x_range = x_max - x_min
                z_range = z_max - z_min
                x_min -= x_range * padding
                x_max += x_range * padding
                z_min -= z_range * padding
                z_max += z_range * padding
                
                x_grid = np.linspace(x_min, x_max, resolution)
                z_grid = np.linspace(z_min, z_max, resolution)
                xx, zz = np.meshgrid(x_grid, z_grid)
                
                mask = np.zeros((resolution, resolution), dtype=bool)
                for i in range(resolution):
                    for j in range(resolution):
                        point = Point(xx[i, j], zz[i, j])
                        if hasattr(self.alpha_shapes['xz'], 'contains'):
                            mask[i, j] = self.alpha_shapes['xz'].contains(point)
                        elif hasattr(self.alpha_shapes['xz'], 'intersects'):
                            mask[i, j] = self.alpha_shapes['xz'].intersects(point)
                
                self.masks['xz'] = {
                    'mask': mask,
                    'extent': [x_min, x_max, z_min, z_max],
                    'resolution': resolution
                }
            
            # YZ mask
            if 'yz' in self.alpha_shapes and self.alpha_shapes['yz'] is not None:
                yz_points = shape_data[:, [1, 2]]
                y_min, y_max = yz_points[:, 0].min(), yz_points[:, 0].max()
                z_min, z_max = yz_points[:, 1].min(), yz_points[:, 1].max()
                
                padding = 0.1
                y_range = y_max - y_min
                z_range = z_max - z_min
                y_min -= y_range * padding
                y_max += y_range * padding
                z_min -= z_range * padding
                z_max += z_range * padding
                
                y_grid = np.linspace(y_min, y_max, resolution)
                z_grid = np.linspace(z_min, z_max, resolution)
                yy, zz = np.meshgrid(y_grid, z_grid)
                
                mask = np.zeros((resolution, resolution), dtype=bool)
                for i in range(resolution):
                    for j in range(resolution):
                        point = Point(yy[i, j], zz[i, j])
                        if hasattr(self.alpha_shapes['yz'], 'contains'):
                            mask[i, j] = self.alpha_shapes['yz'].contains(point)
                        elif hasattr(self.alpha_shapes['yz'], 'intersects'):
                            mask[i, j] = self.alpha_shapes['yz'].intersects(point)
                
                self.masks['yz'] = {
                    'mask': mask,
                    'extent': [y_min, y_max, z_min, z_max],
                    'resolution': resolution
                }
            
            print("Masks computed successfully")
            
        except Exception as e:
            print(f"Mask computation error: {e}")

    def update_plots(self):
        """Update all projection plots"""
        self.update_xy_plot()
        self.update_xz_plot()
        self.update_yz_plot()
        self.update_volume_plot()

    def update_xy_plot(self):
        """Update XY projection plot"""
        self.xy_ax.clear()
        
        if self.current_spots is not None and self.use_spots_for_shape.get():
            self.xy_ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], 
                              c='blue', s=20, alpha=0.6, label=f'Spots ({len(self.current_spots)})')
        
        if self.current_skeleton is not None and self.use_skeleton_for_shape.get():
            self.xy_ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], 
                              c='red', s=10, alpha=0.8, label=f'Skeleton ({len(self.current_skeleton)})')
        
        # Plot alpha shape
        if 'xy' in self.alpha_shapes and self.alpha_shapes['xy'] is not None:
            try:
                if hasattr(self.alpha_shapes['xy'], 'exterior'):
                    # Polygon
                    x, y = self.alpha_shapes['xy'].exterior.xy
                    self.xy_ax.plot(x, y, 'g-', linewidth=2, label='Alpha Shape')
                    self.xy_ax.fill(x, y, alpha=0.2, color='green')
                elif hasattr(self.alpha_shapes['xy'], 'coords'):
                    # LineString or Point
                    x, y = zip(*self.alpha_shapes['xy'].coords)
                    self.xy_ax.plot(x, y, 'g-', linewidth=2, label='Alpha Shape')
            except Exception as e:
                print(f"Error plotting XY alpha shape: {e}")
        
        # Plot mask
        if 'xy' in self.masks:
            mask_data = self.masks['xy']
            extent = mask_data['extent']
            self.xy_ax.imshow(mask_data['mask'], extent=extent, alpha=0.3, 
                             cmap='Reds', origin='lower')
        
        self.xy_ax.set_xlabel('X (μm)')
        self.xy_ax.set_ylabel('Y (μm)')
        self.xy_ax.set_title(f'XY Projection - {self.current_image} - {self.current_cell}')
        self.xy_ax.legend()
        self.xy_ax.grid(True, alpha=0.3)
        self.xy_ax.axis('equal')
        self.xy_canvas.draw()

    def update_xz_plot(self):
        """Update XZ projection plot"""
        self.xz_ax.clear()
        
        if self.current_spots is not None and self.use_spots_for_shape.get():
            self.xz_ax.scatter(self.current_spots[:, 0], self.current_spots[:, 2], 
                              c='blue', s=20, alpha=0.6, label=f'Spots ({len(self.current_spots)})')
        
        if self.current_skeleton is not None and self.use_skeleton_for_shape.get():
            self.xz_ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 2], 
                              c='red', s=10, alpha=0.8, label=f'Skeleton ({len(self.current_skeleton)})')
        
        # Plot alpha shape
        if 'xz' in self.alpha_shapes and self.alpha_shapes['xz'] is not None:
            try:
                if hasattr(self.alpha_shapes['xz'], 'exterior'):
                    x, z = self.alpha_shapes['xz'].exterior.xy
                    self.xz_ax.plot(x, z, 'g-', linewidth=2, label='Alpha Shape')
                    self.xz_ax.fill(x, z, alpha=0.2, color='green')
                elif hasattr(self.alpha_shapes['xz'], 'coords'):
                    x, z = zip(*self.alpha_shapes['xz'].coords)
                    self.xz_ax.plot(x, z, 'g-', linewidth=2, label='Alpha Shape')
            except Exception as e:
                print(f"Error plotting XZ alpha shape: {e}")
        
        # Plot mask
        if 'xz' in self.masks:
            mask_data = self.masks['xz']
            extent = mask_data['extent']
            self.xz_ax.imshow(mask_data['mask'], extent=extent, alpha=0.3, 
                             cmap='Reds', origin='lower')
        
        self.xz_ax.set_xlabel('X (μm)')
        self.xz_ax.set_ylabel('Z (μm)')
        self.xz_ax.set_title(f'XZ Projection - {self.current_image} - {self.current_cell}')
        self.xz_ax.legend()
        self.xz_ax.grid(True, alpha=0.3)
        self.xz_ax.axis('equal')
        self.xz_canvas.draw()

    def update_yz_plot(self):
        """Update YZ projection plot"""
        self.yz_ax.clear()
        
        if self.current_spots is not None and self.use_spots_for_shape.get():
            self.yz_ax.scatter(self.current_spots[:, 1], self.current_spots[:, 2], 
                              c='blue', s=20, alpha=0.6, label=f'Spots ({len(self.current_spots)})')
        
        if self.current_skeleton is not None and self.use_skeleton_for_shape.get():
            self.yz_ax.scatter(self.current_skeleton[:, 1], self.current_skeleton[:, 2], 
                              c='red', s=10, alpha=0.8, label=f'Skeleton ({len(self.current_skeleton)})')
        
        # Plot alpha shape
        if 'yz' in self.alpha_shapes and self.alpha_shapes['yz'] is not None:
            try:
                if hasattr(self.alpha_shapes['yz'], 'exterior'):
                    y, z = self.alpha_shapes['yz'].exterior.xy
                    self.yz_ax.plot(y, z, 'g-', linewidth=2, label='Alpha Shape')
                    self.yz_ax.fill(y, z, alpha=0.2, color='green')
                elif hasattr(self.alpha_shapes['yz'], 'coords'):
                    y, z = zip(*self.alpha_shapes['yz'].coords)
                    self.yz_ax.plot(y, z, 'g-', linewidth=2, label='Alpha Shape')
            except Exception as e:
                print(f"Error plotting YZ alpha shape: {e}")
        
        # Plot mask
        if 'yz' in self.masks:
            mask_data = self.masks['yz']
            extent = mask_data['extent']
            self.yz_ax.imshow(mask_data['mask'], extent=extent, alpha=0.3, 
                             cmap='Reds', origin='lower')
        
        self.yz_ax.set_xlabel('Y (μm)')
        self.yz_ax.set_ylabel('Z (μm)')
        self.yz_ax.set_title(f'YZ Projection - {self.current_image} - {self.current_cell}')
        self.yz_ax.legend()
        self.yz_ax.grid(True, alpha=0.3)
        self.yz_ax.axis('equal')
        self.yz_canvas.draw()

    def save_alpha_shapes(self):
        """Save alpha shapes to files"""
        if not self.alpha_shapes or not self.current_image or not self.current_cell:
            messagebox.showwarning("Warning", "No alpha shapes to save")
            return
        
        try:
            output_dir = os.path.join(self.base_dir, 'alpha_shape_results')
            os.makedirs(output_dir, exist_ok=True)
            
            alpha_val = self.alpha_value.get()
            filename_base = f'{self.current_image}_{self.current_cell}_alpha_{alpha_val:.2f}'
            
            # Save each projection
            for projection, shape in self.alpha_shapes.items():
                if shape is not None:
                    # Save coordinates
                    if hasattr(shape, 'exterior'):
                        coords = list(shape.exterior.coords)
                    elif hasattr(shape, 'coords'):
                        coords = list(shape.coords)
                    else:
                        continue
                    
                    coords_df = pd.DataFrame(coords, columns=[f'{projection[0]}_coord', f'{projection[1]}_coord'])
                    coords_file = os.path.join(output_dir, f'{filename_base}_{projection}_coords.csv')
                    coords_df.to_csv(coords_file, index=False)
                    
                    print(f"Saved {projection} coordinates: {coords_file}")
            
            messagebox.showinfo("Success", f"Alpha shapes saved to:\n{output_dir}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save alpha shapes: {str(e)}")

    def save_masks(self):
        """Save binary masks to files"""
        if not self.masks or not self.current_image or not self.current_cell:
            messagebox.showwarning("Warning", "No masks to save")
            return
        
        try:
            output_dir = os.path.join(self.base_dir, 'alpha_shape_masks')
            os.makedirs(output_dir, exist_ok=True)
            
            alpha_val = self.alpha_value.get()
            resolution = int(self.resolution.get())
            filename_base = f'{self.current_image}_{self.current_cell}_alpha_{alpha_val:.2f}_res_{resolution}'
            
            for projection, mask_data in self.masks.items():
                mask = mask_data['mask'].astype(np.uint8) * 255
                mask_file = os.path.join(output_dir, f'{filename_base}_{projection}_mask.png')
                cv2.imwrite(mask_file, mask)
                
                # Save mask metadata
                metadata = {
                    'projection': projection,
                    'alpha_value': alpha_val,
                    'resolution': resolution,
                    'extent': mask_data['extent'],
                    'image': self.current_image,
                    'cell': self.current_cell
                }
                metadata_file = os.path.join(output_dir, f'{filename_base}_{projection}_metadata.json')
                with open(metadata_file, 'w') as f:
                    json.dump(metadata, f, indent=2)
                
                print(f"Saved {projection} mask: {mask_file}")
            
            messagebox.showinfo("Success", f"Masks saved to:\n{output_dir}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save masks: {str(e)}")

    def batch_process(self):
        """Process all available cells in batch"""
        alpha_val = simpledialog.askfloat("Alpha Value", "Enter alpha value for batch processing:", 
                                         minvalue=0.01, initialvalue=self.alpha_value.get())
        if alpha_val is None:
            return
        
        self.alpha_value.set(alpha_val)
        
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Batch Processing")
        progress_label = ttk.Label(progress_win, text="Processing...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        results = []
        
        for image_name, cells in self.available_images.items():
            for cell_name in cells.keys():
                try:
                    # Load data
                    self.load_cell_data(image_name, cell_name)
                    
                    # Compute alpha shapes
                    self.compute_alpha_shapes()
                    
                    # Collect results
                    shape_data = self.get_shape_data()
                    if shape_data is not None:
                        result = {
                            'image': image_name,
                            'cell': cell_name,
                            'alpha_value': alpha_val,
                            'num_points': len(shape_data),
                            'use_spots': self.use_spots_for_shape.get(),
                            'use_skeleton': self.use_skeleton_for_shape.get()
                        }
                        
                        # Add alpha shape info
                        for projection in ['xy', 'xz', 'yz']:
                            if projection in self.alpha_shapes and self.alpha_shapes[projection] is not None:
                                try:
                                    if hasattr(self.alpha_shapes[projection], 'area'):
                                        result[f'{projection}_area'] = self.alpha_shapes[projection].area
                                    elif hasattr(self.alpha_shapes[projection], 'length'):
                                        result[f'{projection}_length'] = self.alpha_shapes[projection].length
                                    result[f'{projection}_computed'] = True
                                except:
                                    result[f'{projection}_computed'] = False
                            else:
                                result[f'{projection}_computed'] = False
                        
                        results.append(result)
                    
                except Exception as e:
                    print(f"Error processing {image_name}-{cell_name}: {e}")
                
                done += 1
                progress_var.set(done / total)
                progress_label.config(text=f"Processing {image_name}-{cell_name}...")
                self.root.update()
        
        # Save batch results
        if results:
            output_dir = os.path.join(self.base_dir, 'alpha_shape_results')
            os.makedirs(output_dir, exist_ok=True)
            
            results_df = pd.DataFrame(results)
            results_file = os.path.join(output_dir, f'batch_alpha_shape_analysis_alpha_{alpha_val:.2f}.csv')
            results_df.to_csv(results_file, index=False)
            
            print(f"Batch results saved: {results_file}")
        
        progress_win.destroy()
        messagebox.showinfo("Batch Processing Complete", 
                           f"Processed {len(results)} cells.\nResults saved to alpha_shape_results folder.")

    def compute_3d_volume(self):
        """Compute 3D volume from three orthogonal alpha shapes using voxelization"""
        if not self.alpha_shapes or len(self.alpha_shapes) < 3:
            messagebox.showwarning("Warning", "Need all three alpha shapes (XY, XZ, YZ) to compute 3D volume")
            return
        
        # Check if all required shapes exist
        required_projections = ['xy', 'xz', 'yz']
        for proj in required_projections:
            if proj not in self.alpha_shapes or self.alpha_shapes[proj] is None:
                messagebox.showwarning("Warning", f"Missing {proj.upper()} alpha shape")
                return
        
        try:
            shape_data = self.get_shape_data()
            if shape_data is None:
                messagebox.showwarning("Warning", "No data available")
                return
            
            # Get bounding box for all data
            x_min, x_max = shape_data[:, 0].min(), shape_data[:, 0].max()
            y_min, y_max = shape_data[:, 1].min(), shape_data[:, 1].max()
            z_min, z_max = shape_data[:, 2].min(), shape_data[:, 2].max()
            
            # Add padding
            padding = 0.1
            x_range = x_max - x_min
            y_range = y_max - y_min
            z_range = z_max - z_min
            x_min -= x_range * padding
            x_max += x_range * padding
            y_min -= y_range * padding
            y_max += y_range * padding
            z_min -= z_range * padding
            z_max += z_range * padding
            
            # Create 3D voxel grid
            voxel_res = int(self.voxel_resolution.get())
            x_grid = np.linspace(x_min, x_max, voxel_res)
            y_grid = np.linspace(y_min, y_max, voxel_res)
            z_grid = np.linspace(z_min, z_max, voxel_res)
            
            voxel_size = ((x_max - x_min) * (y_max - y_min) * (z_max - z_min)) / (voxel_res ** 3)
            
            print(f"Computing 3D volume with {voxel_res}³ voxels...")
            print(f"Voxel size: {voxel_size:.6f} μm³")
            
            # Initialize volume mask
            volume_mask = np.ones((voxel_res, voxel_res, voxel_res), dtype=bool)
            
            # Process each projection
            projections = {
                'xy': (x_grid, y_grid, 2),  # z-axis index
                'xz': (x_grid, z_grid, 1),  # y-axis index  
                'yz': (y_grid, z_grid, 0)   # x-axis index
            }
            
            for proj_name, (grid1, grid2, axis_idx) in projections.items():
                alpha_shape = self.alpha_shapes[proj_name]
                
                if hasattr(alpha_shape, 'contains') or hasattr(alpha_shape, 'intersects'):
                    # Create 2D mask for this projection
                    proj_mask = np.zeros((len(grid1), len(grid2)), dtype=bool)
                    
                    for i, coord1 in enumerate(grid1):
                        for j, coord2 in enumerate(grid2):
                            point = Point(coord1, coord2)
                            try:
                                if hasattr(alpha_shape, 'contains'):
                                    inside = alpha_shape.contains(point)
                                else:
                                    inside = alpha_shape.intersects(point)
                                proj_mask[i, j] = inside
                            except:
                                proj_mask[i, j] = False
                    
                    # Extrude 2D mask to 3D along the missing axis
                    if axis_idx == 0:  # YZ projection, extrude along X
                        for x_idx in range(voxel_res):
                            volume_mask[x_idx, :, :] &= proj_mask
                    elif axis_idx == 1:  # XZ projection, extrude along Y
                        for y_idx in range(voxel_res):
                            volume_mask[:, y_idx, :] &= proj_mask
                    elif axis_idx == 2:  # XY projection, extrude along Z
                        for z_idx in range(voxel_res):
                            volume_mask[:, :, z_idx] &= proj_mask
            
            # Store results
            self.volume_3d = volume_mask
            volume_voxels = np.sum(volume_mask)
            self.volume_value = volume_voxels * voxel_size
            
            print(f"3D Volume computed successfully:")
            print(f"  Volume voxels: {volume_voxels}")
            print(f"  Total volume: {self.volume_value:.6f} μm³")
            
            # Update display
            self.update_volume_plot()
            
            # Show results
            messagebox.showinfo("3D Volume Computed", 
                              f"Volume: {self.volume_value:.6f} μm³\n"
                              f"Voxels: {volume_voxels}/{voxel_res**3}\n"
                              f"Voxel size: {voxel_size:.6f} μm³")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to compute 3D volume: {str(e)}")
            print(f"3D volume computation error: {e}")
            import traceback
            traceback.print_exc()

    def create_3d_polyhedron(self):
        """Create a 3D polyhedron directly from point cloud data using convex hull or 3D alpha shape"""
        try:
            from scipy.spatial import ConvexHull
            
            # Get the data used for shape generation
            shape_data = self.get_shape_data()
            if shape_data is None or len(shape_data) < 4:
                return None, None
            
            # For very large datasets, sample points to improve performance
            if len(shape_data) > 500:
                sample_idx = np.random.choice(len(shape_data), 500, replace=False)
                shape_data = shape_data[sample_idx]
            
            if ALPHASHAPE_AVAILABLE:
                try:
                    # Try to use 3D alpha shape if alphashape supports 3D
                    alpha_value = self.alpha_value.get()
                    
                    # Check if we can create 3D alpha shape
                    alpha_shape_3d = alphashape.alphashape(shape_data, alpha_value)
                    
                    # Check if we get a valid 3D alpha shape
                    if hasattr(alpha_shape_3d, 'exterior') and alpha_shape_3d.exterior is not None:
                        # For 3D alpha shapes, fall back to convex hull for face generation
                        raise ValueError("3D alpha shape boundary extraction not implemented, using convex hull")
                    else:
                        raise ValueError("3D alpha shape failed, using convex hull")
                        
                except Exception as e:
                    print(f"3D alpha shape failed: {e}, falling back to convex hull")
                    
                    # Fallback to 3D convex hull
                    hull_3d = ConvexHull(shape_data)
                    return hull_3d.points, hull_3d.simplices
            else:
                # Use convex hull directly if alphashape is not available
                print("Using convex hull for 3D polyhedron (alphashape not available)")
                hull_3d = ConvexHull(shape_data)
                return hull_3d.points, hull_3d.simplices
                
        except Exception as e:
            print(f"Error creating 3D polyhedron: {e}")
            return None, None

    def update_volume_plot(self):
        """Update 3D volume visualization with smooth polyhedron boundary"""
        self.volume_ax.clear()
        
        # Use direct 3D polyhedron approach
        vertices, faces = self.create_3d_polyhedron()
        
        if vertices is not None and faces is not None and len(faces) > 0:
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            
            # Validate face indices
            max_vertex_idx = len(vertices) - 1
            valid_faces = []
            for face in faces:
                if all(idx <= max_vertex_idx for idx in face):
                    valid_faces.append(face)
            
            if valid_faces:
                # Create faces for the polyhedron
                mesh = Poly3DCollection(vertices[valid_faces], alpha=0.4, facecolor='lightblue', 
                                      edgecolor='blue', linewidth=0.5)
                self.volume_ax.add_collection3d(mesh)
                surface_label = f'3D Polyhedron ({len(valid_faces)} faces)'
            else:
                surface_label = 'Polyhedron has no valid faces'
            
            # Set equal aspect ratio
            if len(vertices) > 0:
                max_range = max(vertices.max(axis=0) - vertices.min(axis=0))
                mid_x = (vertices[:, 0].max() + vertices[:, 0].min()) * 0.5
                mid_y = (vertices[:, 1].max() + vertices[:, 1].min()) * 0.5
                mid_z = (vertices[:, 2].max() + vertices[:, 2].min()) * 0.5
                
                self.volume_ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
                self.volume_ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
                self.volume_ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)
                
                # Calculate approximate volume using convex hull volume
                try:
                    from scipy.spatial import ConvexHull
                    hull = ConvexHull(vertices)
                    self.volume_value = hull.volume
                except:
                    self.volume_value = 0
            else:
                surface_label = 'No vertices generated'
                
        else:
            surface_label = 'Failed to create polyhedron'
        
        # Also plot original data for reference  
        if self.current_spots is not None and self.use_spots_for_shape.get():
            # Use original coordinates
            spots_x = self.current_spots[:, 0]
            spots_y = self.current_spots[:, 1]
            spots_z = self.current_spots[:, 2]
            
            # Sample spots if too many
            if len(spots_x) > 150:
                sample_idx = np.random.choice(len(spots_x), 150, replace=False)
                spots_x = spots_x[sample_idx]
                spots_y = spots_y[sample_idx] 
                spots_z = spots_z[sample_idx]
            
            self.volume_ax.scatter(spots_x, spots_y, spots_z, 
                                 c='blue', s=8, alpha=0.4, label='Original Spots (sample)')
        
        # Always plot skeleton for reference (regardless of whether it's used for shape)
        if self.current_skeleton is not None:
            # Use original coordinates
            skel_x = self.current_skeleton[:, 0]
            skel_y = self.current_skeleton[:, 1]
            skel_z = self.current_skeleton[:, 2]
            
            # Sample skeleton if too many points
            if len(skel_x) > 100:
                sample_idx = np.random.choice(len(skel_x), 100, replace=False)
                skel_x = skel_x[sample_idx]
                skel_y = skel_y[sample_idx]
                skel_z = skel_z[sample_idx]
            
            self.volume_ax.scatter(skel_x, skel_y, skel_z, 
                                 c='green', s=6, alpha=0.7, label='Skeleton (sample)')
        
        # Set axis labels in micrometers
        self.volume_ax.set_xlabel('X (μm)')
        self.volume_ax.set_ylabel('Y (μm)')
        self.volume_ax.set_zlabel('Z (μm)')
        
        title = f'3D Polyhedron - {self.current_image} - {self.current_cell}'
        if self.volume_value > 0:
            title += f'\nVolume: {self.volume_value:.6f} μm³'
        self.volume_ax.set_title(title)
        
        if 'surface_label' in locals():
            # Create a manual legend entry for the surface
            from matplotlib.lines import Line2D
            legend_elements = []
            if 'surface_label' in locals():
                legend_elements.append(Line2D([0], [0], marker='s', color='w', 
                                             markerfacecolor='lightblue', markersize=10, alpha=0.4,
                                             label=surface_label))
            if self.current_spots is not None and self.use_spots_for_shape.get():
                legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                             markerfacecolor='blue', markersize=8, alpha=0.4,
                                             label='Original Spots'))
            if self.current_skeleton is not None:
                legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                             markerfacecolor='green', markersize=6, alpha=0.7,
                                             label='Skeleton'))
            if legend_elements:
                self.volume_ax.legend(handles=legend_elements)
        
        self.volume_canvas.draw()

    def save_3d_volume(self):
        """Save 3D volume data and visualization"""
        if self.volume_3d is None or not self.current_image or not self.current_cell:
            messagebox.showwarning("Warning", "No 3D volume to save")
            return
        
        try:
            output_dir = os.path.join(self.base_dir, 'alpha_shape_volumes')
            os.makedirs(output_dir, exist_ok=True)
            
            alpha_val = self.alpha_value.get()
            voxel_res = int(self.voxel_resolution.get())
            filename_base = f'{self.current_image}_{self.current_cell}_alpha_{alpha_val:.2f}_voxel_{voxel_res}'
            
            # Save volume data as numpy array
            volume_file = os.path.join(output_dir, f'{filename_base}_volume.npy')
            np.save(volume_file, self.volume_3d)
            
            # Save volume coordinates
            filled_indices = np.where(self.volume_3d)
            if len(filled_indices[0]) > 0:
                volume_coords = np.column_stack([filled_indices[0], filled_indices[1], filled_indices[2]])
                coords_df = pd.DataFrame(volume_coords, columns=['x_voxel', 'y_voxel', 'z_voxel'])
                coords_file = os.path.join(output_dir, f'{filename_base}_coordinates.csv')
                coords_df.to_csv(coords_file, index=False)
            
            # Save metadata
            metadata = {
                'image': self.current_image,
                'cell': self.current_cell,
                'alpha_value': alpha_val,
                'voxel_resolution': voxel_res,
                'volume_cubic_micrometers': float(self.volume_value),
                'volume_voxels': int(np.sum(self.volume_3d)),
                'total_voxels': int(self.volume_3d.size),
                'use_spots': self.use_spots_for_shape.get(),
                'use_skeleton': self.use_skeleton_for_shape.get()
            }
            metadata_file = os.path.join(output_dir, f'{filename_base}_metadata.json')
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
            
            # Save visualization
            viz_file = os.path.join(output_dir, f'{filename_base}_visualization.png')
            self.volume_fig.savefig(viz_file, dpi=300, bbox_inches='tight')
            
            print(f"Saved 3D volume: {volume_file}")
            print(f"Saved coordinates: {coords_file}")
            print(f"Saved metadata: {metadata_file}")
            
            messagebox.showinfo("Success", f"3D volume saved to:\n{output_dir}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save 3D volume: {str(e)}")

    def batch_volume_analysis(self):
        """Batch analysis of 3D volumes for all cells"""
        alpha_val = simpledialog.askfloat("Alpha Value", "Enter alpha value for batch volume analysis:", 
                                         minvalue=0.01, initialvalue=self.alpha_value.get())
        if alpha_val is None:
            return
        
        voxel_res = simpledialog.askinteger("Voxel Resolution", "Enter voxel resolution for batch analysis:", 
                                           minvalue=20, maxvalue=100, initialvalue=self.voxel_resolution.get())
        if voxel_res is None:
            return
        
        self.alpha_value.set(alpha_val)
        self.voxel_resolution.set(voxel_res)
        
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Batch Volume Analysis")
        progress_label = ttk.Label(progress_win, text="Computing volumes...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        results = []
        all_volumes = []
        
        for image_name, cells in self.available_images.items():
            for cell_name in cells.keys():
                try:
                    # Load data
                    self.load_cell_data(image_name, cell_name)
                    
                    # Compute alpha shapes
                    self.compute_alpha_shapes()
                    
                    # Compute 3D volume
                    if len(self.alpha_shapes) >= 3:
                        self.compute_3d_volume()
                        
                        if self.volume_value > 0:
                            result = {
                                'image': image_name,
                                'cell': cell_name,
                                'alpha_value': alpha_val,
                                'voxel_resolution': voxel_res,
                                'volume_cubic_micrometers': self.volume_value,
                                'volume_voxels': int(np.sum(self.volume_3d)) if self.volume_3d is not None else 0,
                                'total_voxels': int(self.volume_3d.size) if self.volume_3d is not None else 0,
                                'volume_fraction': np.sum(self.volume_3d) / self.volume_3d.size if self.volume_3d is not None else 0,
                                'use_spots': self.use_spots_for_shape.get(),
                                'use_skeleton': self.use_skeleton_for_shape.get()
                            }
                            
                            # Add alpha shape areas for reference
                            for projection in ['xy', 'xz', 'yz']:
                                if projection in self.alpha_shapes and self.alpha_shapes[projection] is not None:
                                    try:
                                        if hasattr(self.alpha_shapes[projection], 'area'):
                                            result[f'{projection}_area'] = self.alpha_shapes[projection].area
                                    except:
                                        result[f'{projection}_area'] = 0
                                else:
                                    result[f'{projection}_area'] = 0
                            
                            results.append(result)
                            all_volumes.append(self.volume_value)
                    
                except Exception as e:
                    print(f"Error processing {image_name}-{cell_name}: {e}")
                
                done += 1
                progress_var.set(done / total)
                progress_label.config(text=f"Processing {image_name}-{cell_name}...")
                self.root.update()
        
        progress_win.destroy()
        
        if not results:
            messagebox.showwarning("Warning", "No volumes computed successfully")
            return
        
        # Save batch results
        output_dir = os.path.join(self.base_dir, 'alpha_shape_volumes')
        os.makedirs(output_dir, exist_ok=True)
        
        results_df = pd.DataFrame(results)
        results_file = os.path.join(output_dir, f'batch_volume_analysis_alpha_{alpha_val:.2f}_voxel_{voxel_res}.csv')
        results_df.to_csv(results_file, index=False)
        
        # Create volume distribution plots
        if all_volumes:
            # Histogram
            plt.figure(figsize=(12, 5))
            
            plt.subplot(1, 2, 1)
            n, bins, patches = plt.hist(all_volumes, bins=20, color='skyblue', edgecolor='black', alpha=0.7)
            plt.xlabel('Volume (μm³)')
            plt.ylabel('Number of cells')
            plt.title(f'3D Volume Distribution\n(α={alpha_val}, voxel_res={voxel_res})')
            plt.axvline(np.mean(all_volumes), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(all_volumes):.6f} μm³')
            plt.axvline(np.median(all_volumes), color='orange', linestyle='--', 
                       label=f'Median: {np.median(all_volumes):.6f} μm³')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # Box plot
            plt.subplot(1, 2, 2)
            plt.boxplot(all_volumes, vert=True, patch_artist=True, 
                       boxprops=dict(facecolor='lightgreen'))
            plt.ylabel('Volume (μm³)')
            plt.title(f'Volume Distribution Boxplot\n(N={len(all_volumes)} cells)')
            
            # Add statistics text
            stats_text = f"""Statistics:
Mean: {np.mean(all_volumes):.6f} μm³
Median: {np.median(all_volumes):.6f} μm³
Std: {np.std(all_volumes):.6f} μm³
Min: {np.min(all_volumes):.6f} μm³
Max: {np.max(all_volumes):.6f} μm³
Cells: {len(all_volumes)}"""
            
            plt.figtext(0.02, 0.5, stats_text, fontsize=10, 
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            plt.tight_layout()
            
            # Save plots
            plot_file = os.path.join(output_dir, f'volume_distribution_alpha_{alpha_val:.2f}_voxel_{voxel_res}.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.show()
            
            print(f"Batch volume analysis completed:")
            print(f"  Results saved: {results_file}")
            print(f"  Plots saved: {plot_file}")
            print(f"  Volume statistics:")
            print(f"    Mean: {np.mean(all_volumes):.6f} μm³")
            print(f"    Median: {np.median(all_volumes):.6f} μm³")
            print(f"    Range: {np.min(all_volumes):.6f} - {np.max(all_volumes):.6f} μm³")
        
        messagebox.showinfo("Batch Volume Analysis Complete", 
                           f"Processed {len(results)} cells.\n"
                           f"Mean volume: {np.mean(all_volumes):.6f} μm³\n"
                           f"Results saved to alpha_shape_volumes folder.")


def main():
    # Check if required packages are installed
    try:
        import alphashape
        import shapely
        import cv2
    except ImportError as e:
        print(f"Missing required package: {e}")
        print("Please install required packages:")
        print("pip install alphashape shapely opencv-python")
        return
    
    root = tk.Tk()
    app = AlphaShapeProjectionAnalyzer(root)
    root.mainloop()

if __name__ == "__main__":
    main() 