import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import threading
import tkinter.simpledialog as simpledialog
import json

class Interactive3DBatchVisualizer:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive 3D Cell Batch Visualization")
        self.root.geometry("1500x950")

        # Path settings
        self.base_dir = r'Y333 ATP6 ATP2'
        self.skeleton_root = os.path.join(self.base_dir, 'extracted_cells')
        # self.spots_root = self.base_dir
        self.channel = 'atp6'
        self.spots_root = os.path.join(self.base_dir, f'{self.channel}_spots')

        # Data structures
        self.available_images = {}  # {image_name: {cell_name: {spots_data, skeleton_file}}}
        self.coordinate_mappings = {} # To store coordinate mappings from json
        self.selected_image = tk.StringVar()
        self.selected_cell = tk.StringVar()
        self.current_cell = None
        self.current_image = None

        # Other variables
        self.original_spots = None
        self.original_skeleton = None
        self.current_skeleton = None
        self.current_spots = None
        self.distances = None
        self.current_skeleton_translation = np.array([0, 0, 0])  # Current skeleton translation
        self.rotation_x = tk.DoubleVar(value=0.0)
        self.rotation_y = tk.DoubleVar(value=0.0)
        self.rotation_z = tk.DoubleVar(value=0.0)
        self.use_y_flip = tk.BooleanVar(value=True)
        self.use_z_flip = tk.BooleanVar(value=False)
        self.auto_compare_yflip = tk.BooleanVar(value=True)
        self.auto_translate_skeleton = tk.BooleanVar(value=False)

        # Create interface
        self.create_interface()
        # Scan all images and cells
        self.scan_available_images_and_cells()
        # Auto load first image and cell
        if self.available_images:
            first_image = list(self.available_images.keys())[0]
            self.selected_image.set(first_image)
            self.on_image_changed()

    def scan_available_images_and_cells(self):
        """Scan all images and cells, spot files are matched with image folders through index and field of view double matching"""
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

            # Load coordinate mapping for the image
            mapping_file = os.path.join(image_path, 'coordinate_mapping.json')
            if os.path.exists(mapping_file):
                with open(mapping_file, 'r') as f:
                    self.coordinate_mappings[image_folder] = json.load(f)
            else:
                print(f"Warning: coordinate_mapping.json not found in {image_path}")
                continue # Skip if no mapping file

            # Parse index and fov
            parts = image_folder.split('_')
            # Assume the last two are index and fov respectively
            if len(parts) < 2:
                continue
            index = parts[-2]
            fov = parts[-1]
            # Match spot file: contains _index_ and _fov_, and has _spots
            matched_spot = None
            for spot_file in spot_files:
                if f'_{index}_' in spot_file and f'_{fov}_' in spot_file and '_spots' in spot_file:
                    matched_spot = os.path.join(self.spots_root, spot_file)
                    break
            if not matched_spot:
                print(f"Spot file not found for {image_folder} (index={index}, fov={fov})")
                continue
            # Parse spot file
            cells_data = self.parse_spots_file(matched_spot)
            # Scan skeleton files
            skeleton_files = [f for f in os.listdir(image_path) if f.startswith('cell_') and f.endswith('.txt')]
            image_cells = {}
            for skel_file in skeleton_files:
                try:
                    cell_num = int(skel_file.split('_')[1].split('.')[0])
                    corresponding_cell_name = None
                    for cell_name in cells_data.keys():
                        if cell_name.startswith('Cell_') and int(cell_name.split('_')[1]) == cell_num:
                            corresponding_cell_name = cell_name
                            break
                    if corresponding_cell_name:
                        # Check if this cell has valid spots data
                        spots_data = cells_data[corresponding_cell_name]['spots']
                        if len(spots_data) > 0:  # Only add cells with spots
                            image_cells[corresponding_cell_name] = {
                                'spots_data': cells_data[corresponding_cell_name],
                                'skeleton_file': os.path.join(image_path, skel_file)
                            }
                        else:
                            print(f"Skipping {image_folder}-{corresponding_cell_name}: no spots data")
                except Exception as e:
                    print(f"Failed to parse skeleton filename {skel_file}: {e}")
            if image_cells:
                self.available_images[image_folder] = image_cells
        # Update image dropdown
        if hasattr(self, 'image_combo'):
            self.image_combo['values'] = list(self.available_images.keys())
            if self.available_images:
                self.selected_image.set(list(self.available_images.keys())[0])

    def parse_spots_file(self, spots_file):
        """Parse spots file and extract cell outlines and spots data"""
        with open(spots_file, 'r') as f:
            lines = f.readlines()
        
        cells_data = {}
        current_cell = None
        reading_spots = False
        reading_x_pos = False
        reading_y_pos = False
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            if line.startswith('CELL'):
                current_cell = line.split('\t')[1]
                cells_data[current_cell] = {'spots': [], 'outline_x': [], 'outline_y': []}
                reading_spots = False
                reading_x_pos = False
                reading_y_pos = False
                continue
            elif line.startswith('X_POS'):
                reading_x_pos = True
                reading_spots = False
                reading_y_pos = False
                # Parse X coordinates from the same line
                parts = line.split('\t')[1:]  # Skip 'X_POS' and get the coordinates
                x_coords = []
                for part in parts:
                    part = part.strip()
                    if part and part != 'END':
                        try:
                            x_coords.append(float(part))
                        except ValueError:
                            pass
                if x_coords:
                    cells_data[current_cell]['outline_x'].extend(x_coords)
                continue
            elif line.startswith('Y_POS'):
                reading_y_pos = True
                reading_spots = False
                reading_x_pos = False
                # Parse Y coordinates from the same line
                parts = line.split('\t')[1:]  # Skip 'Y_POS' and get the coordinates
                y_coords = []
                for part in parts:
                    part = part.strip()
                    if part and part != 'END':
                        try:
                            y_coords.append(float(part))
                        except ValueError:
                            pass
                if y_coords:
                    cells_data[current_cell]['outline_y'].extend(y_coords)
                continue
            elif line == 'SPOTS':
                reading_spots = True
                reading_x_pos = False
                reading_y_pos = False
                continue
            elif line == 'END':
                reading_spots = False
                reading_x_pos = False
                reading_y_pos = False
                continue
            
            # Parse X_POS coordinates - now handle multiple numbers in one line
            if reading_x_pos and line and current_cell and not line.startswith('CELL'):
                try:
                    parts = line.split('\t')
                    x_coords = []
                    for part in parts:
                        part = part.strip()
                        if part and part != 'END':
                            try:
                                x_coords.append(float(part))
                            except ValueError:
                                pass
                    if x_coords:
                        cells_data[current_cell]['outline_x'].extend(x_coords)
                except Exception as e:
                    print(f"Error parsing X_POS line: {line[:50]}... Error: {e}")
                    continue
            
            # Parse Y_POS coordinates - now handle multiple numbers in one line
            elif reading_y_pos and line and current_cell and not line.startswith('CELL'):
                try:
                    parts = line.split('\t')
                    y_coords = []
                    for part in parts:
                        part = part.strip()
                        if part and part != 'END':
                            try:
                                y_coords.append(float(part))
                            except ValueError:
                                pass
                    if y_coords:
                        cells_data[current_cell]['outline_y'].extend(y_coords)
                except Exception as e:
                    print(f"Error parsing Y_POS line: {line[:50]}... Error: {e}")
                    continue
            
            # Parse SPOTS data
            elif reading_spots and line and not line.startswith('Pos_Y') and current_cell:
                if line.startswith('CELL'):
                    current_cell = line.split('\t')[1]
                    cells_data[current_cell] = {'spots': [], 'outline_x': [], 'outline_y': []}
                    reading_spots = False
                else:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        try:
                            y = float(parts[0])
                            x = float(parts[1])
                            z = float(parts[2])
                            cells_data[current_cell]['spots'].append([z, y, x])
                        except ValueError:
                            continue
        

        
        return cells_data

    def read_skeleton_data(self, skeleton_file, mapping_data=None, cell_name=None, pixel_size_xy=0.0645):
        """
        完全按照analyze_alignment.py中load_skeleton_txt方法加载skeleton数据并应用坐标转换
        """
        df = pd.read_csv(skeleton_file, sep='\t')
        # 提取坐标列 (x, y, z) - 这些是相对于截取区域的坐标
        coords = df[['x', 'y', 'z']].values
        print(f"Skeleton数据加载完成: {len(coords)} 个点")
        print(f"相对坐标范围: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
        
        # 如果提供了mapping数据，转换为绝对坐标（完全复制analyze_alignment.py中的逻辑）
        if mapping_data and cell_name:
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                x_offset = crop_info['x_offset']
                y_offset = crop_info['y_offset']
                
                # 将像素偏移量转换为微米单位，并加到相对坐标上
                offset_x_um = x_offset * pixel_size_xy
                offset_y_um = y_offset * pixel_size_xy
                
                coords[:, 0] += offset_x_um  # X坐标加偏移
                coords[:, 1] += offset_y_um  # Y坐标加偏移
                # Z坐标不需要偏移，因为是3D图像的深度方向
                
                print(f"应用坐标偏移: X+{offset_x_um:.3f}μm (像素{x_offset}), Y+{offset_y_um:.3f}μm (像素{y_offset})")
                print(f"绝对坐标范围: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
            else:
                print(f"警告: 在mapping文件中未找到{cell_name}")
        else:
            if mapping_data:
                print(f"警告: mapping数据存在但cell_name为空")
            else:
                print("警告: 未提供mapping数据")
            print("使用相对坐标（未应用偏移）")
        
        return coords

    def get_cell_outline_data(self, image_name, cell_name):
        """Get cell outline data from parsed spots file"""
        if image_name not in self.available_images or cell_name not in self.available_images[image_name]:
            return None
        
        try:
            cell_info = self.available_images[image_name][cell_name]['spots_data']
            outline_x = cell_info.get('outline_x', [])
            outline_y = cell_info.get('outline_y', [])
            
            if len(outline_x) == 0 or len(outline_y) == 0:
                return None
            
            # Make sure x and y coordinates have the same length
            min_len = min(len(outline_x), len(outline_y))
            if min_len == 0:
                return None
                
            # Create 2D outline coordinates (x, y, z=0 for now)
            outline_coords = np.zeros((min_len, 3))
            outline_coords[:, 0] = outline_x[:min_len]
            outline_coords[:, 1] = outline_y[:min_len]
            outline_coords[:, 2] = 0  # Set Z to 0 for 2D outline
            
            return outline_coords
            
        except Exception as e:
            print(f"Error getting outline data for {image_name}-{cell_name}: {e}")
            return None

    def convert_and_scale_coordinates(self, spots_px, skeleton_coords, mapping_data, pixel_size_xy=64.5, pixel_size_z=200.0):
        """
        完全按照analyze_alignment.py中load_spots_fishquant方法的正确逻辑进行坐标转换
        
        关键：
        1. spots_px来自parse_spots_file，格式为[z, y, x]（像素单位）
        2. 但analyze_alignment.py中从文件读取的是[y, x, z]格式
        3. 需要转换为纳米→Y轴翻转→重排为(x,y,z)→转微米
        """
        if spots_px.ndim != 2 or spots_px.shape[1] != 3:
            raise ValueError(f"Invalid spots data format: expected 2D array with 3 columns, got shape {spots_px.shape}")

        # 第一步：将spots_px[z,y,x]转换为analyze_alignment.py中期望的[y,x,z]格式
        # spots_px格式：[z, y, x] -> coords格式：[y, x, z]
        coords = np.zeros_like(spots_px)
        coords[:, 0] = spots_px[:, 1]  # Y
        coords[:, 1] = spots_px[:, 2]  # X  
        coords[:, 2] = spots_px[:, 0]  # Z
        
        # 第二步：转换为纳米单位（完全复制analyze_alignment.py的逻辑）
        coords_nm = coords.copy()
        coords_nm[:, 0] *= pixel_size_xy  # Y轴: 像素 -> 纳米
        coords_nm[:, 1] *= pixel_size_xy  # X轴: 像素 -> 纳米  
        coords_nm[:, 2] *= pixel_size_z   # Z轴: 像素 -> 纳米

        # 第三步：Y轴翻转（在纳米单位上进行，完全按照analyze_alignment.py的方法）
        if self.use_y_flip.get() and mapping_data:
            # 使用mapping数据中的crop区域信息进行精确翻转
            crop_info = mapping_data['crop_region']
            y_start = crop_info['y_start']  # 像素坐标
            y_end = crop_info['y_end']      # 像素坐标
            
            # 转换为纳米坐标 (像素 * 像素大小)
            y_start_nm = y_start * pixel_size_xy
            y_end_nm = y_end * pixel_size_xy
            
            # 基于crop区域进行翻转：new_y = (y_start + y_end) - old_y
            flip_center = y_start_nm + y_end_nm
            coords_nm[:, 0] = flip_center - coords_nm[:, 0]  # 在纳米单位上翻转Y坐标
            
            print(f"Y轴翻转: 基于crop区域 y_start={y_start_nm:.1f}nm, y_end={y_end_nm:.1f}nm")
        elif self.use_y_flip.get():
            print("警告: 需要mapping数据进行精确Y轴翻转，使用简单翻转")
            if len(coords_nm) > 0:
                y_min, y_max = coords_nm[:, 0].min(), coords_nm[:, 0].max()
                flip_center = y_min + y_max
                coords_nm[:, 0] = flip_center - coords_nm[:, 0]

        # 第四步：重新排列为(x, y, z)格式以匹配skeleton（完全复制analyze_alignment.py）
        # coords_nm格式: [y, x, z] -> spots_nm_xyz格式: [x, y, z]
        spots_nm_xyz = coords_nm[:, [1, 0, 2]]  # from (y, x, z) to (x, y, z)

        # 第五步：转换为微米单位（就像analyze_alignment.py: spots_coords_um = spots_coords / 1000.0）
        spots_um_xyz = spots_nm_xyz / 1000.0

        print(f"Spots坐标转换完成: X=[{spots_um_xyz[:,0].min():.3f}, {spots_um_xyz[:,0].max():.3f}], Y=[{spots_um_xyz[:,1].min():.3f}, {spots_um_xyz[:,1].max():.3f}], Z=[{spots_um_xyz[:,2].min():.3f}, {spots_um_xyz[:,2].max():.3f}]")

        # 第六步：应用Z轴翻转（如果需要）
        if self.use_z_flip.get():
            z_center = (spots_um_xyz[:, 2].max() + spots_um_xyz[:, 2].min()) / 2
            spots_um_xyz[:, 2] = 2 * z_center - spots_um_xyz[:, 2]

        return {
            'original': spots_um_xyz,
            'flip_y': spots_um_xyz, # For compatibility
            'scale_x': 1.0,
            'scale_y': 1.0
        }

    def create_interface(self):
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        title_label = ttk.Label(control_frame, text="Batch Cell Analysis Control", font=('Arial', 14, 'bold'))
        title_label.pack(pady=(0, 20))
        self.create_image_selection(control_frame)
        self.create_cell_selection(control_frame)
        self.create_transform_options(control_frame)
        self.create_rotation_controls(control_frame)
        self.create_action_buttons(control_frame)
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.create_3d_plot(plot_frame)
        self.create_result_area(plot_frame)

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
        transform_frame = ttk.LabelFrame(parent, text="Coordinate Transform Options", padding=10)
        transform_frame.pack(fill=tk.X, pady=(0, 15))
        y_flip_check = ttk.Checkbutton(transform_frame, text="Y-axis Flip", variable=self.use_y_flip, command=self.on_transform_change)
        y_flip_check.pack(anchor=tk.W)
        z_flip_check = ttk.Checkbutton(transform_frame, text="Z-axis Flip (along XY plane)", variable=self.use_z_flip, command=self.on_transform_change)
        z_flip_check.pack(anchor=tk.W)
        
        # Add separator
        ttk.Separator(transform_frame, orient='horizontal').pack(fill=tk.X, pady=(10, 10))
        
        # Batch analysis options
        ttk.Label(transform_frame, text="Batch Analysis Options:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        auto_compare_check = ttk.Checkbutton(transform_frame, text="Y-flip", variable=self.auto_compare_yflip)
        auto_compare_check.pack(anchor=tk.W)
        auto_translate_check = ttk.Checkbutton(transform_frame, text="Use aligned skeleton", variable=self.auto_translate_skeleton)
        auto_translate_check.pack(anchor=tk.W)
        
        ttk.Label(transform_frame, text="(Z flip: around cell range center)", font=('Arial', 8), foreground='gray').pack(anchor=tk.W)

    def create_rotation_controls(self, parent):
        rotation_frame = ttk.LabelFrame(parent, text="Skeleton Rotation Control", padding=10)
        rotation_frame.pack(fill=tk.X, pady=(0, 15))
        x_frame = ttk.Frame(rotation_frame)
        x_frame.pack(fill=tk.X, pady=2)
        ttk.Label(x_frame, text="X-axis:", width=6).pack(side=tk.LEFT)
        x_scale = ttk.Scale(x_frame, from_=-180, to=180, variable=self.rotation_x, orient=tk.HORIZONTAL, command=self.on_rotation_change)
        x_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        x_value = ttk.Label(x_frame, text="0°", width=6)
        x_value.pack(side=tk.RIGHT)
        self.x_value_label = x_value
        y_frame = ttk.Frame(rotation_frame)
        y_frame.pack(fill=tk.X, pady=2)
        ttk.Label(y_frame, text="Y-axis:", width=6).pack(side=tk.LEFT)
        y_scale = ttk.Scale(y_frame, from_=-180, to=180, variable=self.rotation_y, orient=tk.HORIZONTAL, command=self.on_rotation_change)
        y_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        y_value = ttk.Label(y_frame, text="0°", width=6)
        y_value.pack(side=tk.RIGHT)
        self.y_value_label = y_value
        z_frame = ttk.Frame(rotation_frame)
        z_frame.pack(fill=tk.X, pady=2)
        ttk.Label(z_frame, text="Z-axis:", width=6).pack(side=tk.LEFT)
        z_scale = ttk.Scale(z_frame, from_=-180, to=180, variable=self.rotation_z, orient=tk.HORIZONTAL, command=self.on_rotation_change)
        z_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        z_value = ttk.Label(z_frame, text="0°", width=6)
        z_value.pack(side=tk.RIGHT)
        self.z_value_label = z_value

    def create_action_buttons(self, parent):
        button_frame = ttk.LabelFrame(parent, text="Actions", padding=10)
        button_frame.pack(fill=tk.X, pady=(0, 15))
        reset_btn = ttk.Button(button_frame, text="Reset Rotation", command=self.reset_rotation)
        reset_btn.pack(fill=tk.X, pady=2)
        calculate_btn = ttk.Button(button_frame, text="Calculate Distance Analysis", command=self.calculate_and_show_results)
        calculate_btn.pack(fill=tk.X, pady=5)
        save_btn = ttk.Button(button_frame, text="Save Results", command=self.save_results)
        save_btn.pack(fill=tk.X, pady=2)
        batch_btn = ttk.Button(button_frame, text="Batch Distance Analysis", command=self.batch_distance_analysis)
        batch_btn.pack(fill=tk.X, pady=5)
        
        outlier_btn = ttk.Button(button_frame, text="Analyse Outliers", command=self.analyse_outliers)
        outlier_btn.pack(fill=tk.X, pady=5)

    def create_3d_plot(self, parent):
        self.fig = Figure(figsize=(10, 8), dpi=100)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def create_result_area(self, parent):
        pass

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
            self.load_cell_data(image_name, cell_name)
            self.update_cell_info()
            self.update_3d_plot()

    def on_transform_change(self):
        if self.current_image and self.current_cell:
            self.load_cell_data(self.current_image, self.current_cell)
            # Reset translation info because coordinate transformation changed
            self.current_skeleton_translation = np.array([0, 0, 0])
            self.update_3d_plot()

    def update_cell_info(self):
        if self.current_cell and hasattr(self, 'cell_info_label'):
            spots_count = len(self.current_spots) if self.current_spots is not None else 0
            skeleton_count = len(self.current_skeleton) if self.current_skeleton is not None else 0
            info_text = f"Spots: {spots_count}, Skeleton: {skeleton_count}"
            self.cell_info_label.config(text=info_text)

    def on_rotation_change(self, value=None):
        self.x_value_label.config(text=f"{self.rotation_x.get():.0f}°")
        self.y_value_label.config(text=f"{self.rotation_y.get():.0f}°")
        self.z_value_label.config(text=f"{self.rotation_z.get():.0f}°")
        self.update_skeleton_rotation()
        self.update_3d_plot()

    def update_skeleton_rotation(self):
        if self.original_skeleton is None:
            return
        rx = np.radians(self.rotation_x.get())
        ry = np.radians(self.rotation_y.get())
        rz = np.radians(self.rotation_z.get())
        rotation = R.from_euler('xyz', [rx, ry, rz])
        center = np.mean(self.original_skeleton, axis=0)
        centered_skeleton = self.original_skeleton - center
        rotated_skeleton = rotation.apply(centered_skeleton)
        self.current_skeleton = rotated_skeleton + center

    def update_3d_plot(self):
        self.ax.clear()
        if self.current_spots is not None and len(self.current_spots) > 0:
            self.ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], self.current_spots[:, 2], c='blue', s=30, alpha=0.6, label=f'Spots ({len(self.current_spots)})')
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            self.ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], self.current_skeleton[:, 2], c='red', s=10, alpha=0.8, label=f'Skeleton ({len(self.current_skeleton)})')
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
        if (self.current_spots is not None and len(self.current_spots) > 0) or (self.current_skeleton is not None and len(self.current_skeleton) > 0):
            self.ax.legend()
        self.set_equal_aspect_3d()
        self.canvas.draw()
        self.update_cell_info()

    def set_equal_aspect_3d(self):
        if self.current_spots is not None and self.current_skeleton is not None:
            all_points = np.vstack([self.current_spots, self.current_skeleton])
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

    def reset_rotation(self):
        self.rotation_x.set(0)
        self.rotation_y.set(0)
        self.rotation_z.set(0)
        self.on_rotation_change()

    def load_cell_data(self, image_name, cell_name):
        """Load specified cell data using the perfect alignment scheme from analyze_alignment.py"""
        if image_name not in self.available_images or cell_name not in self.available_images[image_name]:
            print(f"Cell {cell_name} in image {image_name} not available")
            return
        
        try:
            self.current_image = image_name
            self.current_cell = cell_name
            cell_info = self.available_images[image_name][cell_name]
            
            # Get coordinate mapping for the current cell (从analyze_alignment.py方案)
            cell_number = int(cell_name.split("_")[1])
            cell_id_str = f"{cell_number:03d}"
            mapping_cell_name = f'cell_{cell_id_str}'
            mapping_data = self.coordinate_mappings.get(image_name, {}).get(mapping_cell_name)
            
            # 获取像素大小（模拟analyze_alignment.py中从spots文件读取的过程）
            pixel_size_xy = 64.5  # nm, 默认值
            pixel_size_z = 200.0  # nm, 默认值
            
            # Load skeleton data with coordinate transformation (使用analyze_alignment.py的方法)
            skeleton_coords = self.read_skeleton_data(
                cell_info['skeleton_file'], 
                mapping_data={mapping_cell_name: mapping_data} if mapping_data else None,
                cell_name=mapping_cell_name,
                pixel_size_xy=pixel_size_xy/1000  # 转换为微米
            )
            
            print(f"Successfully loaded skeleton data for {cell_name}:")
            print(f"  Skeleton coordinate range: X=[{skeleton_coords[:,0].min():.3f}, {skeleton_coords[:,0].max():.3f}], Y=[{skeleton_coords[:,1].min():.3f}, {skeleton_coords[:,1].max():.3f}], Z=[{skeleton_coords[:,2].min():.3f}, {skeleton_coords[:,2].max():.3f}]")
            
            self.original_skeleton = skeleton_coords.copy()
            self.current_skeleton = skeleton_coords.copy()
            
            # Get spots coordinates (in pixels, format: z, y, x)
            spots_px = np.array(cell_info['spots_data']['spots'])
            
            # Check spots data format
            if spots_px.ndim != 2 or spots_px.shape[1] != 3 or spots_px.shape[0] == 0:
                print(f"Warning: Invalid or empty spots data for {image_name}-{cell_name}: shape {spots_px.shape}")
                self.current_spots = None
                self.original_spots = None
            else:
                # Transform coordinates using the verified alignment scheme
                conversion_results = self.convert_and_scale_coordinates(
                    spots_px, skeleton_coords, mapping_data, 
                    pixel_size_xy=pixel_size_xy, pixel_size_z=pixel_size_z
                )
                
                spots_um_xyz = conversion_results['original'].copy()
                
                # 关键的中心校正步骤（就像analyze_alignment.py中的最后步骤）
                skeleton_center = np.mean(skeleton_coords, axis=0)
                spots_center = np.mean(spots_um_xyz, axis=0)
                center_diff = skeleton_center - spots_center
                center_distance = np.linalg.norm(center_diff)
                
                print(f"中心校正分析:")
                print(f"  Skeleton中心: X={skeleton_center[0]:.3f}, Y={skeleton_center[1]:.3f}, Z={skeleton_center[2]:.3f}")
                print(f"  Spots中心: X={spots_center[0]:.3f}, Y={spots_center[1]:.3f}, Z={spots_center[2]:.3f}")
                print(f"  中心距离: {center_distance:.3f}μm")
                
                # 应用中心校正（analyze_alignment.py: spots_corrected = spots_coords_um - center_diff）
                # 注意：analyze_alignment.py中是减法，因为是spots向skeleton靠近
                spots_corrected = spots_um_xyz - center_diff
                
                self.current_spots = spots_corrected.copy()
                self.original_spots = spots_corrected.copy()
                
                print(f"校正后spots数据:")
                print(f"  Spots coordinate range: X=[{self.current_spots[:,0].min():.3f}, {self.current_spots[:,0].max():.3f}], Y=[{self.current_spots[:,1].min():.3f}, {self.current_spots[:,1].max():.3f}], Z=[{self.current_spots[:,2].min():.3f}, {self.current_spots[:,2].max():.3f}]")
            
            # Process outline data using same transformation parameters
            self.current_outline = self.process_outline_data(image_name, cell_name, skeleton_coords)
            
            # Update window title
            self.root.title(f"Interactive 3D Cell Batch Visualization - {image_name} - {cell_name}")
            
            spots_count = len(self.original_spots) if self.original_spots is not None else 0
            outline_count = len(self.current_outline) if self.current_outline is not None else 0
            print(f"Successfully loaded {image_name}-{cell_name} data:")
            print(f"  Spots count: {spots_count}")
            print(f"  Skeleton points count: {len(self.original_skeleton)}")
            print(f"  Outline points count: {outline_count}")
            
            # Reset rotation
            self.reset_rotation()
            
        except Exception as e:
            print(f"Failed to load {image_name}-{cell_name} data: {e}")
            import traceback
            traceback.print_exc()
    
    def process_outline_data(self, image_name, cell_name, skeleton_coords):
        """
        Transforms raw outline coordinates (in pixels, relative to the original image)
        to the skeleton's micrometer-based coordinate system.
        """
        voxel_size_y = 0.160  # um/pixel
        voxel_size_x = 0.160  # um/pixel

        # Get raw outline coordinates (x, y) in pixels
        outline_px = self.get_cell_outline_data(image_name, cell_name)
        if outline_px is None or len(outline_px) == 0:
            return None

        # Get coordinate mapping for the current cell
        mapping_data = self.coordinate_mappings.get(image_name, {}).get(f'cell_{int(cell_name.split("_")[1]):03d}')

        # Apply pixel offset
        if mapping_data:
            x_offset = mapping_data['crop_region']['x_offset']
            y_offset = mapping_data['crop_region']['y_offset']
            outline_px[:, 0] -= x_offset
            outline_px[:, 1] -= y_offset

        # Convert pixel coordinates to micrometers
        outline_um = outline_px.copy()
        outline_um[:, 0] = outline_um[:, 0] * voxel_size_x
        outline_um[:, 1] = outline_um[:, 1] * voxel_size_y

        # Apply Y-axis flip (if enabled), which should be consistent with spots
        if self.use_y_flip.get():
            outline_um[:, 1] = -outline_um[:, 1]

        return outline_um

    def calculate_distances(self):
        if self.current_spots is None or self.current_skeleton is None:
            return np.array([])
        if len(self.current_spots) == 0 or len(self.current_skeleton) == 0:
            return np.array([])
        distances = cdist(self.current_spots, self.current_skeleton)
        min_distances = np.min(distances, axis=1)
        return min_distances

    def optimize_skeleton_translation(self, spots, skeleton, search_range=2.0, step_size=0.2):
        """
        Automatically find optimal skeleton translation position to minimize median distance to spots
        
        Args:
            spots: spots coordinates (N, 3)
            skeleton: skeleton coordinates (M, 3)  
            search_range: search range in each direction (μm)
            step_size: search step size (μm)
        
        Returns:
            dict: {'optimal_skeleton': optimal translated skeleton, 'translation': optimal translation vector, 'median_distance': minimum median distance}
        """
        if len(spots) == 0 or len(skeleton) == 0:
            return {'optimal_skeleton': skeleton, 'translation': np.array([0, 0, 0]), 'median_distance': float('inf')}
        
        from scipy.optimize import minimize
        
        def objective_function(translation):
            """Objective function: return median distance from translated skeleton to spots"""
            translated_skeleton = skeleton + translation.reshape(1, 3)
            distances = cdist(spots, translated_skeleton)
            min_distances = np.min(distances, axis=1)
            return np.median(min_distances)
        
        # Initial guess: zero translation
        initial_translation = np.array([0.0, 0.0, 0.0])
        
        # Set search bounds
        bounds = [(-search_range, search_range) for _ in range(3)]
        
        # Use L-BFGS-B algorithm for optimization
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
            print(f"Optimization failed, using original skeleton: {e}")
            distances = cdist(spots, skeleton)
            min_distances = np.min(distances, axis=1)
            median_distance = np.median(min_distances)
            return {
                'optimal_skeleton': skeleton,
                'translation': np.array([0, 0, 0]),
                'median_distance': median_distance
            }

    def calculate_and_show_results(self):
        if self.current_spots is None or self.current_skeleton is None:
            messagebox.showwarning("Warning", "Please select a cell first")
            return
        
        # If auto translation is enabled, perform translation optimization first
        if self.auto_translate_skeleton.get():
            translation_result = self.optimize_skeleton_translation(self.current_spots, self.current_skeleton)
            self.current_skeleton = translation_result['optimal_skeleton']
            self.current_skeleton_translation = translation_result['translation']
            print(f"Auto translation optimization completed, translation: X={self.current_skeleton_translation[0]:.3f}, Y={self.current_skeleton_translation[1]:.3f}, Z={self.current_skeleton_translation[2]:.3f}")
        else:
            self.current_skeleton_translation = np.array([0, 0, 0])
        
        self.distances = self.calculate_distances()
        if len(self.distances) == 0:
            messagebox.showwarning("Warning", "Cannot calculate distances, please check data")
            return
        self.create_result_window()
        print(f"{self.current_image} - {self.current_cell} calculation completed:")
        print(f"  Mean distance: {np.mean(self.distances):.2f} μm")
        print(f"  Median distance: {np.median(self.distances):.2f} μm")
        print(f"  Standard deviation: {np.std(self.distances):.2f} μm")

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
        
        # Use processed outline data
        cell_outline = self.current_outline if hasattr(self, 'current_outline') else None
        
        # Draw cell outline first (so it appears behind other elements)
        if cell_outline is not None and len(cell_outline) > 0:
            # Plot outline as connected lines (close the polygon)
            outline_2d = cell_outline[:, :2]  # Use only X, Y coordinates
            # Close the polygon by adding the first point at the end
            closed_outline = np.vstack([outline_2d, outline_2d[0]])
            ax.plot(closed_outline[:, 0], closed_outline[:, 1], 'k-', linewidth=2, alpha=0.8, label='Cell outline')
        
        # Plot spots and skeleton as before
        if self.current_spots is not None and len(self.current_spots) > 0:
            scatter = ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], c=self.distances, cmap='viridis', s=50, alpha=0.7)
            fig.colorbar(scatter, ax=ax, label='Distance to skeleton (μm)')
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], c='red', s=10, alpha=0.5, label='Skeleton')
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
        if (self.current_skeleton is not None and len(self.current_skeleton) > 0) or cell_outline is not None:
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
            ax.axvline(np.mean(self.distances), color='red', linestyle='--', label=f'Mean: {np.mean(self.distances):.2f} μm')
            ax.axvline(np.median(self.distances), color='orange', linestyle='--', label=f'Median: {np.median(self.distances):.2f} μm')
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
        ax.grid(True, alpha=0.3)
        
        if len(self.distances) > 0:
            transform_status = []
            if self.use_y_flip.get():
                transform_status.append('Y-flip')
            if self.use_z_flip.get():
                transform_status.append('Z-flip')
            if self.auto_translate_skeleton.get():
                transform_status.append('Auto-translate')
            transform_text = '+'.join(transform_status) if transform_status else 'Original'
            
            translation_text = ""
            if self.auto_translate_skeleton.get():
                translation_text = f"\nTranslation: X={self.current_skeleton_translation[0]:.3f} Y={self.current_skeleton_translation[1]:.3f} Z={self.current_skeleton_translation[2]:.3f}"
            
            stats_text = f"""Statistics:\nCount: {len(self.distances)}\nMean: {np.mean(self.distances):.3f} μm\nMedian: {np.median(self.distances):.3f} μm\nStd: {np.std(self.distances):.3f} μm\nMin: {np.min(self.distances):.3f} μm\nMax: {np.max(self.distances):.3f} μm\n\nTransform: {transform_text}\nRotation: X={self.rotation_x.get():.0f}° Y={self.rotation_y.get():.0f}° Z={self.rotation_z.get():.0f}°{translation_text}"""
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)
        canvas = FigureCanvasTkAgg(fig, parent_frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas.draw()

    def save_results(self):
        if self.distances is None or len(self.distances) == 0:
            messagebox.showwarning("Warning", "No results to save, please calculate distance analysis first")
            return
        if not self.current_cell or not self.current_image:
            messagebox.showwarning("Warning", "No cell or image selected")
            return
        try:
            output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
            os.makedirs(output_dir, exist_ok=True)
            detail_data = pd.DataFrame({
                'spot_x_um': self.current_spots[:, 0],
                'spot_y_um': self.current_spots[:, 1],
                'spot_z_um': self.current_spots[:, 2],
                'distance_to_skeleton_um': self.distances,
                'rotation_x': self.rotation_x.get(),
                'rotation_y': self.rotation_y.get(),
                'rotation_z': self.rotation_z.get(),
                'y_flip_used': self.use_y_flip.get(),
                'z_flip_used': self.use_z_flip.get(),
                'auto_compare_enabled': self.auto_compare_yflip.get(),
                'auto_translate_enabled': self.auto_translate_skeleton.get(),
                'skeleton_translation_x': self.current_skeleton_translation[0],
                'skeleton_translation_y': self.current_skeleton_translation[1],
                'skeleton_translation_z': self.current_skeleton_translation[2],
                'image_name': self.current_image,
                'cell_name': self.current_cell
            })
            detail_file = os.path.join(output_dir, f'{self.current_image}_{self.current_cell}_interactive_analysis.csv')
            detail_data.to_csv(detail_file, index=False)
            skeleton_data = pd.DataFrame({
                'skeleton_x_um': self.current_skeleton[:, 0],
                'skeleton_y_um': self.current_skeleton[:, 1],
                'skeleton_z_um': self.current_skeleton[:, 2],
            })
            skeleton_file = os.path.join(output_dir, f'{self.current_image}_{self.current_cell}_rotated_skeleton.csv')
            skeleton_data.to_csv(skeleton_file, index=False)
            summary_data = {
                'image_name': self.current_image,
                'cell_name': self.current_cell,
                'num_spots': len(self.current_spots),
                'num_skeleton_points': len(self.current_skeleton),
                'y_flip_used': self.use_y_flip.get(),
                'z_flip_used': self.use_z_flip.get(),
                'auto_compare_enabled': self.auto_compare_yflip.get(),
                'auto_translate_enabled': self.auto_translate_skeleton.get(),
                'skeleton_translation_x': self.current_skeleton_translation[0],
                'skeleton_translation_y': self.current_skeleton_translation[1],
                'skeleton_translation_z': self.current_skeleton_translation[2],
                'rotation_x': self.rotation_x.get(),
                'rotation_y': self.rotation_y.get(),
                'rotation_z': self.rotation_z.get(),
                'mean_distance': np.mean(self.distances),
                'median_distance': np.median(self.distances),
                'std_distance': np.std(self.distances),
                'min_distance': np.min(self.distances),
                'max_distance': np.max(self.distances)
            }
            summary_df = pd.DataFrame([summary_data])
            summary_file = os.path.join(output_dir, f'{self.current_image}_{self.current_cell}_summary.csv')
            summary_df.to_csv(summary_file, index=False)
            messagebox.showinfo("Success", f"Results saved to:\n{output_dir}")
        except Exception as e:
            messagebox.showerror("Error", f"Save failed: {str(e)}")

    def batch_distance_analysis(self):
        """Batch analysis of distance distribution for all images and all cells, output aggregated histogram and boxplot, and calculate threshold ratio statistics, with progress bar"""
        threshold = simpledialog.askfloat("Threshold Input", "Please enter distance threshold (μm):", minvalue=0.0, initialvalue=0.6)
        if threshold is None:
            return
        # Create progress bar window
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Batch Analysis Progress")
        progress_label = ttk.Label(progress_win, text="Analyzing...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        all_distances = []
        summary_rows = []
        cell_spot_counts = []
        all_spots_details = []  # Record detailed information for all spots
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        for image_name, cells in self.available_images.items():
            for cell_name, cell_info in cells.items():
                spots_nm = np.array(cell_info['spots_data']['spots'])
                
                # Get coordinate mapping for the current cell (使用analyze_alignment.py方案)
                cell_number = int(cell_name.split("_")[1])
                cell_id_str = f"{cell_number:03d}"
                mapping_cell_name = f'cell_{cell_id_str}'
                mapping_data = self.coordinate_mappings.get(image_name, {}).get(mapping_cell_name)
                
                # 加载skeleton数据，应用coordinate transformation
                pixel_size_xy = 64.5  # nm
                skeleton_coords = self.read_skeleton_data(
                    cell_info['skeleton_file'], 
                    mapping_data={mapping_cell_name: mapping_data} if mapping_data else None,
                    cell_name=mapping_cell_name,
                    pixel_size_xy=pixel_size_xy/1000  # 转换为微米
                )
                # Skip cells without spots or skeleton (should not happen due to filtering during scanning)
                if spots_nm.ndim != 2 or spots_nm.shape[0] == 0 or skeleton_coords.shape[0] == 0:
                    done += 1
                    progress_var.set(done / total)
                    self.root.update()
                    continue
                try:
                    # 使用统一的坐标转换方法，完全复制load_cell_data中的逻辑
                    if self.auto_compare_yflip.get():
                        # Auto comparison mode: calculate y_flip True/False separately, choose optimal
                        best_result = None
                        best_mean = None
                        original_y_flip = self.use_y_flip.get()  # 保存原始设置
                        
                        for y_flip in [True, False]:
                            self.use_y_flip.set(y_flip)
                            # 使用完全相同的坐标转换逻辑
                            conversion_results = self.convert_and_scale_coordinates(
                                spots_nm, skeleton_coords, mapping_data, 
                                pixel_size_xy=pixel_size_xy, pixel_size_z=200.0
                            )
                            spots_um_xyz = conversion_results['original']
                            
                            # 应用中心校正（完全复制analyze_alignment.py和load_cell_data中的逻辑）
                            skeleton_center = np.mean(skeleton_coords, axis=0)
                            spots_center = np.mean(spots_um_xyz, axis=0)
                            center_diff = skeleton_center - spots_center
                            spots_corrected = spots_um_xyz - center_diff
                            
                            spots = spots_corrected
                            skeleton = skeleton_coords
                            
                            # Apply auto translation optimization (if enabled)
                            skeleton_translation = np.array([0, 0, 0])
                            if self.auto_translate_skeleton.get():
                                translation_result = self.optimize_skeleton_translation(spots, skeleton)
                                skeleton = translation_result['optimal_skeleton']
                                skeleton_translation = translation_result['translation']
                            
                            distances = cdist(spots, skeleton)
                            min_distances = np.min(distances, axis=1)
                            mean_dist = np.mean(min_distances)
                            
                            if (best_mean is None) or (mean_dist < best_mean):
                                best_mean = mean_dist
                                best_result = {
                                    'spots': spots.copy(),
                                    'skeleton': skeleton.copy(),
                                    'min_distances': min_distances.copy(),
                                    'y_flip': y_flip,
                                    'skeleton_translation': skeleton_translation.copy()
                                }
                        
                        # Use optimal result
                        spots = best_result['spots']
                        skeleton = best_result['skeleton']
                        min_distances = best_result['min_distances']
                        y_flip_used = best_result['y_flip']
                        skeleton_translation = best_result['skeleton_translation']
                        
                        # 恢复原始设置
                        self.use_y_flip.set(original_y_flip)
                        
                    else:
                        # Use current interface settings, no comparison
                        conversion_results = self.convert_and_scale_coordinates(
                            spots_nm, skeleton_coords, mapping_data, 
                            pixel_size_xy=pixel_size_xy, pixel_size_z=200.0
                        )
                        spots_um_xyz = conversion_results['original']
                        
                        # 应用中心校正（完全复制analyze_alignment.py和load_cell_data中的逻辑）
                        skeleton_center = np.mean(skeleton_coords, axis=0)
                        spots_center = np.mean(spots_um_xyz, axis=0)
                        center_diff = skeleton_center - spots_center
                        spots_corrected = spots_um_xyz - center_diff
                        
                        spots = spots_corrected
                        skeleton = skeleton_coords
                        y_flip_used = self.use_y_flip.get()
                        
                        # Apply auto translation optimization (if enabled)
                        skeleton_translation = np.array([0, 0, 0])  # Default no translation
                        if self.auto_translate_skeleton.get():
                            translation_result = self.optimize_skeleton_translation(spots, skeleton)
                            skeleton = translation_result['optimal_skeleton']
                            skeleton_translation = translation_result['translation']
                        
                        # Calculate distances
                        distances = cdist(spots, skeleton)
                        min_distances = np.min(distances, axis=1)
                    below = (min_distances < threshold)
                    for i, d in enumerate(min_distances):
                        all_distances.append(d)
                        # Record detailed information for each spot
                        all_spots_details.append({
                            'image': image_name,
                            'cell': cell_name,
                            'spot_index': i,
                            'spot_x_um': spots[i, 0],
                            'spot_y_um': spots[i, 1], 
                            'spot_z_um': spots[i, 2],
                            'distance_to_skeleton_um': d,
                            'exceeds_threshold': d >= threshold,
                            'threshold_um': threshold,
                            'y_flip_used': y_flip_used if self.auto_compare_yflip.get() else self.use_y_flip.get(),
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
                    print(f"Error in {image_name}-{cell_name}: {e}")
                done += 1
                progress_var.set(done / total)
                self.root.update()
        # Save all distance distribution
        output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
        os.makedirs(output_dir, exist_ok=True)
        import pandas as pd
        dist_df = pd.DataFrame({'distance': all_distances})
        dist_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_distances.csv'), index=False)
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_summary.csv'), index=False)
        
        # Save all spots' detailed information
        if all_spots_details:
            spots_details_df = pd.DataFrame(all_spots_details)
            spots_details_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_spots_details.csv'), index=False)
            
            # Save spots exceeding threshold
            exceeds_threshold_df = spots_details_df[spots_details_df['exceeds_threshold'] == True]
            exceeds_threshold_df.to_csv(os.path.join(output_dir, f'{self.channel}_spots_exceed_threshold_{threshold}.csv'), index=False)
            
            print(f"Total analysis: {len(spots_details_df)} spots")
            print(f"Spots exceeding threshold {threshold} μm: {len(exceeds_threshold_df)} ({len(exceeds_threshold_df)/len(spots_details_df)*100:.1f}%)")
        # Plot aggregated histogram and boxplot
        if not dist_df.empty:
            import matplotlib.pyplot as plt
            mean_val = dist_df['distance'].mean()
            median_val = dist_df['distance'].median()
            below_ratio = (dist_df['distance'] < threshold).mean()
            cell_count = len(cell_spot_counts)
            mean_spots_per_cell = np.mean(cell_spot_counts) if cell_spot_counts else 0
            plt.figure(figsize=(10, 6))
            n, bins, patches = plt.hist(dist_df['distance'], bins=30, color='skyblue', edgecolor='black', alpha=0.7)
            plt.xlabel('Distance to skeleton (μm)')
            plt.ylabel('Number of molecules')
            plt.title(f'{self.channel.upper()} - All Cells Distance Distribution Histogram (N={len(dist_df)})')
            # Plot mean, median, threshold vertical lines
            plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
            plt.axvline(median_val, color='orange', linestyle='--', label=f'Median: {median_val:.2f}')
            plt.axvline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            # Annotate ratio below threshold, total cells, average spots
            plt.text(0.98, 0.95, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().transAxes, ha='right', va='top', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_histogram.png'))
            plt.close()

            plt.figure(figsize=(8, 6))
            plt.boxplot(dist_df['distance'], vert=True, patch_artist=True, boxprops=dict(facecolor='lightgreen'))
            plt.ylabel('Distance to skeleton (μm)')
            plt.title(f'{self.channel.upper()} - All Cells Distance Boxplot (N={len(dist_df)})')
            # Plot threshold line
            plt.axhline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            # Annotate ratio below threshold, total cells, average spots
            plt.text(1.05, threshold, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().get_yaxis_transform(which='grid'), ha='left', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_boxplot.png'))
            plt.close()
        progress_win.destroy()
        exceed_count = len(exceeds_threshold_df) if 'exceeds_threshold_df' in locals() else 0
        total_spots = len(all_spots_details) if all_spots_details else 0
        messagebox.showinfo("Batch Analysis Completed", f"Distance distribution histogram and boxplot for all cells saved to:\n{output_dir}\n\nTotal spots: {total_spots}\nSpots exceeding threshold {threshold} μm: {exceed_count}")

    def analyse_outliers(self):
        """分析outliers并在新窗口中显示包含outliers的图像和细胞"""
        # 询问阈值
        threshold = simpledialog.askfloat("Outlier阈值输入", "请输入outlier距离阈值 (μm)：", minvalue=0.0, initialvalue=0.6)
        if threshold is None:
            return
        
        # 创建进度条窗口
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Analyse Outliers Progress")
        progress_label = ttk.Label(progress_win, text="Analyzing outliers...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        
        # 查找有outliers的图像和细胞
        outlier_images = {}  # {image_name: {cell_name: outlier_info}}
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        
        for image_name, cells in self.available_images.items():
            for cell_name, cell_info in cells.items():
                spots_nm = np.array(cell_info['spots_data']['spots'])
                
                # Get coordinate mapping for the current cell (使用analyze_alignment.py方案)
                cell_number = int(cell_name.split("_")[1])
                cell_id_str = f"{cell_number:03d}"
                mapping_cell_name = f'cell_{cell_id_str}'
                mapping_data = self.coordinate_mappings.get(image_name, {}).get(mapping_cell_name)
                
                # 加载skeleton数据，应用coordinate transformation
                pixel_size_xy = 64.5  # nm
                skeleton_coords = self.read_skeleton_data(
                    cell_info['skeleton_file'], 
                    mapping_data={mapping_cell_name: mapping_data} if mapping_data else None,
                    cell_name=mapping_cell_name,
                    pixel_size_xy=pixel_size_xy/1000  # 转换为微米
                )
                
                # 跳过无spot或无skeleton的cell（由于扫描时已经过滤，这里应该不会有无spots的情况）
                if spots_nm.ndim != 2 or spots_nm.shape[0] == 0 or skeleton_coords.shape[0] == 0:
                    done += 1
                    progress_var.set(done / total)
                    self.root.update()
                    continue
                
                try:
                    # 使用完全相同的坐标转换逻辑（复制load_cell_data的方法）
                    conversion_results = self.convert_and_scale_coordinates(
                        spots_nm, skeleton_coords, mapping_data, 
                        pixel_size_xy=pixel_size_xy, pixel_size_z=200.0
                    )
                    spots_um_xyz = conversion_results['original']
                    
                    # 应用中心校正（完全复制analyze_alignment.py和load_cell_data中的逻辑）
                    skeleton_center = np.mean(skeleton_coords, axis=0)
                    spots_center = np.mean(spots_um_xyz, axis=0)
                    center_diff = skeleton_center - spots_center
                    spots_corrected = spots_um_xyz - center_diff
                    
                    spots = spots_corrected
                    skeleton = skeleton_coords
                    
                    # 应用自动平移优化（如果启用）
                    skeleton_translation = np.array([0, 0, 0])
                    if self.auto_translate_skeleton.get():
                        translation_result = self.optimize_skeleton_translation(spots, skeleton)
                        skeleton = translation_result['optimal_skeleton']
                        skeleton_translation = translation_result['translation']
                    
                    distances = cdist(spots, skeleton)
                    min_distances = np.min(distances, axis=1)
                    
                    # 检查是否有outliers
                    outlier_indices = np.where(min_distances >= threshold)[0]
                    if len(outlier_indices) > 0:
                        outlier_images[image_name] = outlier_images.get(image_name, {})
                        outlier_images[image_name][cell_name] = {
                            'spots': spots,
                            'skeleton': skeleton,
                            'distances': min_distances,
                            'outlier_indices': outlier_indices,
                            'threshold': threshold,
                            'skeleton_translation': skeleton_translation,
                            'original_cell_info': cell_info
                        }
                
                except Exception as e:
                    print(f"Error analyzing {image_name}-{cell_name}: {e}")
                
                done += 1
                progress_var.set(done / total)
                self.root.update()
        
        progress_win.destroy()
        
        if not outlier_images:
            messagebox.showinfo("分析结果", f"没有找到超出阈值 {threshold} μm 的outliers")
            return
        
        # 统计信息
        total_outlier_cells = sum(len(cells) for cells in outlier_images.values())
        total_outliers = sum(len(cell_data['outlier_indices']) for image_cells in outlier_images.values() 
                           for cell_data in image_cells.values())
        
        # 创建outlier分析窗口
        self.create_outlier_analysis_window(outlier_images, threshold, total_outlier_cells, total_outliers)

    def create_outlier_analysis_window(self, outlier_images, threshold, total_outlier_cells, total_outliers):
        """创建outlier分析窗口"""
        outlier_window = tk.Toplevel(self.root)
        outlier_window.title(f"Outlier Analysis - Threshold: {threshold} μm")
        outlier_window.geometry("1500x950")
        
        # 创建主框架
        main_frame = ttk.Frame(outlier_window)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 控制面板
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # 标题
        title_label = ttk.Label(control_frame, text=f"Outlier Analysis\nThreshold: {threshold} μm", 
                               font=('Arial', 14, 'bold'))
        title_label.pack(pady=(0, 20))
        
        # 统计信息
        stats_frame = ttk.LabelFrame(control_frame, text="Statistics", padding=10)
        stats_frame.pack(fill=tk.X, pady=(0, 15))
        
        stats_text = f"Images with outliers: {len(outlier_images)}\nCells with outliers: {total_outlier_cells}\nTotal outlier spots: {total_outliers}"
        stats_label = ttk.Label(stats_frame, text=stats_text, font=('Arial', 10))
        stats_label.pack(anchor=tk.W)
        
        # 图像选择
        image_frame = ttk.LabelFrame(control_frame, text="Image Selection", padding=10)
        image_frame.pack(fill=tk.X, pady=(0, 15))
        
        outlier_selected_image = tk.StringVar()
        ttk.Label(image_frame, text="Select Image:").pack(anchor=tk.W)
        image_combo = ttk.Combobox(image_frame, textvariable=outlier_selected_image, 
                                  values=list(outlier_images.keys()), state="readonly", width=20)
        image_combo.pack(fill=tk.X, pady=5)
        
        # 细胞选择
        cell_frame = ttk.LabelFrame(control_frame, text="Cell Selection", padding=10)
        cell_frame.pack(fill=tk.X, pady=(0, 15))
        
        outlier_selected_cell = tk.StringVar()
        ttk.Label(cell_frame, text="Select Cell:").pack(anchor=tk.W)
        cell_combo = ttk.Combobox(cell_frame, textvariable=outlier_selected_cell, 
                                 values=[], state="readonly", width=15)
        cell_combo.pack(fill=tk.X, pady=5)
        
        # 细胞信息标签
        cell_info_label = ttk.Label(cell_frame, text="", font=('Arial', 9))
        cell_info_label.pack(anchor=tk.W, pady=(5, 0))
        
        # 3D图形
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        outlier_fig = Figure(figsize=(10, 8), dpi=100)
        outlier_ax = outlier_fig.add_subplot(111, projection='3d')
        outlier_canvas = FigureCanvasTkAgg(outlier_fig, plot_frame)
        outlier_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        def update_outlier_cell_selection(event=None):
            """更新细胞选择下拉框"""
            image_name = outlier_selected_image.get()
            if image_name and image_name in outlier_images:
                cell_names = list(outlier_images[image_name].keys())
                cell_combo['values'] = cell_names
                if cell_names:
                    outlier_selected_cell.set(cell_names[0])
                    update_outlier_plot()
        
        def update_outlier_plot(event=None):
            """更新outlier 3D图"""
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
            
            # 绘制normal spots (蓝色)
            normal_indices = np.where(distances < threshold_val)[0]
            if len(normal_indices) > 0:
                normal_spots = spots[normal_indices]
                outlier_ax.scatter(normal_spots[:, 0], normal_spots[:, 1], normal_spots[:, 2], 
                                 c='blue', s=30, alpha=0.6, label=f'Normal spots ({len(normal_indices)})')
            
            # 绘制outlier spots (红色)
            if len(outlier_indices) > 0:
                outlier_spots = spots[outlier_indices]
                outlier_ax.scatter(outlier_spots[:, 0], outlier_spots[:, 1], outlier_spots[:, 2], 
                                 c='red', s=50, alpha=0.8, label=f'Outlier spots ({len(outlier_indices)})')
            
            # 绘制skeleton (灰色)
            if len(skeleton) > 0:
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
            
            outlier_ax.set_title(f'{self.channel.upper()} - {image_name} - {cell_name} - Outlier Analysis{flip_status}\nThreshold: {threshold_val} μm')
            outlier_ax.legend()
            
            # 设置等比例
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
            
            # 更新细胞信息
            total_spots = len(spots)
            outlier_count = len(outlier_indices)
            normal_count = len(normal_indices)
            outlier_ratio = outlier_count / total_spots * 100 if total_spots > 0 else 0
            
            info_text = f"Total: {total_spots}, Normal: {normal_count}, Outliers: {outlier_count} ({outlier_ratio:.1f}%)"
            cell_info_label.config(text=info_text)
        
        # 绑定事件
        image_combo.bind('<<ComboboxSelected>>', update_outlier_cell_selection)
        cell_combo.bind('<<ComboboxSelected>>', update_outlier_plot)
        
        # 自动选择第一个图像和细胞
        if outlier_images:
            first_image = list(outlier_images.keys())[0]
            outlier_selected_image.set(first_image)
            update_outlier_cell_selection()
        
        # 旋转控制
        rotation_frame = ttk.LabelFrame(control_frame, text="Skeleton Rotation Control", padding=10)
        rotation_frame.pack(fill=tk.X, pady=(0, 15))
        
        # 旋转变量
        outlier_rotation_x = tk.DoubleVar(value=0.0)
        outlier_rotation_y = tk.DoubleVar(value=0.0)
        outlier_rotation_z = tk.DoubleVar(value=0.0)
        
        # X轴旋转
        x_frame = ttk.Frame(rotation_frame)
        x_frame.pack(fill=tk.X, pady=2)
        ttk.Label(x_frame, text="X-axis:", width=6).pack(side=tk.LEFT)
        x_scale = ttk.Scale(x_frame, from_=-180, to=180, variable=outlier_rotation_x, orient=tk.HORIZONTAL)
        x_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        outlier_x_value_label = ttk.Label(x_frame, text="0°", width=6)
        outlier_x_value_label.pack(side=tk.RIGHT)
        
        # Y轴旋转
        y_frame = ttk.Frame(rotation_frame)
        y_frame.pack(fill=tk.X, pady=2)
        ttk.Label(y_frame, text="Y-axis:", width=6).pack(side=tk.LEFT)
        y_scale = ttk.Scale(y_frame, from_=-180, to=180, variable=outlier_rotation_y, orient=tk.HORIZONTAL)
        y_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        outlier_y_value_label = ttk.Label(y_frame, text="0°", width=6)
        outlier_y_value_label.pack(side=tk.RIGHT)
        
        # Z轴旋转
        z_frame = ttk.Frame(rotation_frame)
        z_frame.pack(fill=tk.X, pady=2)
        ttk.Label(z_frame, text="Z-axis:", width=6).pack(side=tk.LEFT)
        z_scale = ttk.Scale(z_frame, from_=-180, to=180, variable=outlier_rotation_z, orient=tk.HORIZONTAL)
        z_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        outlier_z_value_label = ttk.Label(z_frame, text="0°", width=6)
        outlier_z_value_label.pack(side=tk.RIGHT)
        
        # 操作按钮
        action_frame = ttk.LabelFrame(control_frame, text="Actions", padding=10)
        action_frame.pack(fill=tk.X, pady=(15, 0))
        
        def reset_outlier_rotation():
            """重置旋转"""
            outlier_rotation_x.set(0)
            outlier_rotation_y.set(0)
            outlier_rotation_z.set(0)
            update_outlier_rotation()
        
        def update_outlier_rotation(value=None):
            """更新旋转显示和应用旋转"""
            outlier_x_value_label.config(text=f"{outlier_rotation_x.get():.0f}°")
            outlier_y_value_label.config(text=f"{outlier_rotation_y.get():.0f}°")
            outlier_z_value_label.config(text=f"{outlier_rotation_z.get():.0f}°")
            
            # 应用旋转到当前选中的细胞
            image_name = outlier_selected_image.get()
            cell_name = outlier_selected_cell.get()
            if image_name in outlier_images and cell_name in outlier_images[image_name]:
                cell_data = outlier_images[image_name][cell_name]
                
                # Get coordinate mapping for proper skeleton loading
                cell_number = int(cell_name.split("_")[1])
                cell_id_str = f"{cell_number:03d}"
                mapping_cell_name = f'cell_{cell_id_str}'
                mapping_data = self.coordinate_mappings.get(image_name, {}).get(mapping_cell_name)
                
                # 加载原始skeleton数据（使用完全相同的方法）
                pixel_size_xy = 64.5  # nm
                original_skeleton = self.read_skeleton_data(
                    cell_data['original_cell_info']['skeleton_file'],
                    mapping_data={mapping_cell_name: mapping_data} if mapping_data else None,
                    cell_name=mapping_cell_name,
                    pixel_size_xy=pixel_size_xy/1000  # 转换为微米
                )
                
                # 应用旋转
                rx = np.radians(outlier_rotation_x.get())
                ry = np.radians(outlier_rotation_y.get())
                rz = np.radians(outlier_rotation_z.get())
                rotation = R.from_euler('xyz', [rx, ry, rz])
                center = np.mean(original_skeleton, axis=0)
                centered_skeleton = original_skeleton - center
                rotated_skeleton = rotation.apply(centered_skeleton)
                
                # 应用平移（如果有）
                final_skeleton = rotated_skeleton + center + cell_data['skeleton_translation']
                
                # 更新cell_data
                outlier_images[image_name][cell_name]['skeleton'] = final_skeleton
                
                # 重新计算距离
                spots = cell_data['spots']
                distances = cdist(spots, final_skeleton)
                min_distances = np.min(distances, axis=1)
                outlier_images[image_name][cell_name]['distances'] = min_distances
                
                # 重新计算outlier indices
                outlier_indices = np.where(min_distances >= threshold)[0]
                outlier_images[image_name][cell_name]['outlier_indices'] = outlier_indices
                
                # 更新图形
                update_outlier_plot()
        
        def calculate_outlier_distances():
            """计算当前选中细胞的距离分析"""
            image_name = outlier_selected_image.get()
            cell_name = outlier_selected_cell.get()
            if not image_name or not cell_name or image_name not in outlier_images or cell_name not in outlier_images[image_name]:
                messagebox.showwarning("Warning", "请先选择一个细胞")
                return
            
            cell_data = outlier_images[image_name][cell_name]
            distances = cell_data['distances']
            
            if len(distances) == 0:
                messagebox.showwarning("Warning", "没有距离数据可分析")
                return
            
            # 创建结果窗口
            self.create_outlier_result_window(image_name, cell_name, cell_data, threshold)
        
        def save_outlier_cell_results():
            """保存当前选中细胞的分析结果"""
            image_name = outlier_selected_image.get()
            cell_name = outlier_selected_cell.get()
            if not image_name or not cell_name or image_name not in outlier_images or cell_name not in outlier_images[image_name]:
                messagebox.showwarning("Warning", "请先选择一个细胞")
                return
            
            try:
                cell_data = outlier_images[image_name][cell_name]
                spots = cell_data['spots']
                skeleton = cell_data['skeleton']
                distances = cell_data['distances']
                outlier_indices = cell_data['outlier_indices']
                
                output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
                os.makedirs(output_dir, exist_ok=True)
                
                # 保存详细数据
                detail_data = pd.DataFrame({
                    'spot_x_um': spots[:, 0],
                    'spot_y_um': spots[:, 1],
                    'spot_z_um': spots[:, 2],
                    'distance_to_skeleton_um': distances,
                    'is_outlier': [i in outlier_indices for i in range(len(distances))],
                    'threshold_um': threshold,
                    'rotation_x': outlier_rotation_x.get(),
                    'rotation_y': outlier_rotation_y.get(),
                    'rotation_z': outlier_rotation_z.get(),
                    'y_flip_used': self.use_y_flip.get(),
                    'z_flip_used': self.use_z_flip.get(),
                    'auto_translate_enabled': self.auto_translate_skeleton.get(),
                    'skeleton_translation_x': cell_data['skeleton_translation'][0],
                    'skeleton_translation_y': cell_data['skeleton_translation'][1],
                    'skeleton_translation_z': cell_data['skeleton_translation'][2],
                    'image_name': image_name,
                    'cell_name': cell_name
                })
                detail_file = os.path.join(output_dir, f'{image_name}_{cell_name}_outlier_analysis.csv')
                detail_data.to_csv(detail_file, index=False)
                
                # 保存skeleton数据
                skeleton_data = pd.DataFrame({
                    'skeleton_x_um': skeleton[:, 0],
                    'skeleton_y_um': skeleton[:, 1],
                    'skeleton_z_um': skeleton[:, 2],
                })
                skeleton_file = os.path.join(output_dir, f'{image_name}_{cell_name}_outlier_rotated_skeleton.csv')
                skeleton_data.to_csv(skeleton_file, index=False)
                
                # 保存摘要
                outlier_count = len(outlier_indices)
                total_spots = len(spots)
                summary_data = {
                    'image_name': image_name,
                    'cell_name': cell_name,
                    'num_spots': total_spots,
                    'num_outlier_spots': outlier_count,
                    'outlier_ratio': outlier_count / total_spots if total_spots > 0 else 0,
                    'num_skeleton_points': len(skeleton),
                    'threshold_um': threshold,
                    'y_flip_used': self.use_y_flip.get(),
                    'z_flip_used': self.use_z_flip.get(),
                    'auto_translate_enabled': self.auto_translate_skeleton.get(),
                    'skeleton_translation_x': cell_data['skeleton_translation'][0],
                    'skeleton_translation_y': cell_data['skeleton_translation'][1],
                    'skeleton_translation_z': cell_data['skeleton_translation'][2],
                    'rotation_x': outlier_rotation_x.get(),
                    'rotation_y': outlier_rotation_y.get(),
                    'rotation_z': outlier_rotation_z.get(),
                    'mean_distance': np.mean(distances),
                    'median_distance': np.median(distances),
                    'std_distance': np.std(distances),
                    'min_distance': np.min(distances),
                    'max_distance': np.max(distances)
                }
                summary_df = pd.DataFrame([summary_data])
                summary_file = os.path.join(output_dir, f'{image_name}_{cell_name}_outlier_summary.csv')
                summary_df.to_csv(summary_file, index=False)
                
                messagebox.showinfo("Save Successful", f"Outlier analysis results saved to:\n{output_dir}")
            except Exception as e:
                messagebox.showerror("Error", f"Save failed: {str(e)}")
        
        def save_outlier_details():
            """Save all outlier detailed information to CSV"""
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
                outlier_file = os.path.join(output_dir, f'{self.channel}_outlier_analysis_{threshold}.csv')
                outlier_df.to_csv(outlier_file, index=False)
                messagebox.showinfo("Save Successful", f"Outlier detailed information saved to:\n{outlier_file}")
        
        # 添加按钮
        reset_btn = ttk.Button(action_frame, text="Reset Rotation", command=reset_outlier_rotation)
        reset_btn.pack(fill=tk.X, pady=2)
        
        calculate_btn = ttk.Button(action_frame, text="Calculate Distance Analysis", command=calculate_outlier_distances)
        calculate_btn.pack(fill=tk.X, pady=5)
        
        save_cell_btn = ttk.Button(action_frame, text="Save Current Cell Results", command=save_outlier_cell_results)
        save_cell_btn.pack(fill=tk.X, pady=2)
        
        save_all_btn = ttk.Button(action_frame, text="Save All Outlier Details", command=save_outlier_details)
        save_all_btn.pack(fill=tk.X, pady=2)
        
        # 绑定旋转事件
        x_scale.configure(command=update_outlier_rotation)
        y_scale.configure(command=update_outlier_rotation)
        z_scale.configure(command=update_outlier_rotation)

    def create_outlier_result_window(self, image_name, cell_name, cell_data, threshold):
        """为outlier分析创建结果窗口"""
        result_window = tk.Toplevel(self.root)
        result_window.title(f"Outlier Analysis Results - {image_name} - {cell_name}")
        result_window.geometry("1200x800")
        
        result_notebook = ttk.Notebook(result_window)
        result_notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # XY投影页面
        xy_frame = ttk.Frame(result_notebook)
        result_notebook.add(xy_frame, text="XY Projection")
        
        # 距离分布页面
        dist_frame = ttk.Frame(result_notebook)
        result_notebook.add(dist_frame, text="Distance Distribution")
        
        # 创建XY投影
        self.create_outlier_xy_projection_in_window(xy_frame, image_name, cell_name, cell_data, threshold)
        
        # 创建距离分布
        self.create_outlier_distance_distribution_in_window(dist_frame, image_name, cell_name, cell_data, threshold)
        
        result_notebook.select(0)
        result_window.lift()
        result_window.focus_force()

    def create_outlier_xy_projection_in_window(self, parent_frame, image_name, cell_name, cell_data, threshold):
        """Create outlier analysis XY projection plot"""
        fig = Figure(figsize=(10, 8), dpi=100)
        ax = fig.add_subplot(111)
        
        spots = cell_data['spots']
        distances = cell_data['distances']
        outlier_indices = cell_data['outlier_indices']
        skeleton = cell_data['skeleton']
        
        # Note: Outline display is temporarily disabled in outlier analysis
        # TODO: Implement proper outline transformation for outlier analysis
        cell_outline = None
        
        # Draw normal spots
        normal_indices = np.where(distances < threshold)[0]
        if len(normal_indices) > 0:
            normal_spots = spots[normal_indices]
            normal_distances = distances[normal_indices]
            scatter_normal = ax.scatter(normal_spots[:, 0], normal_spots[:, 1], 
                                      c=normal_distances, cmap='viridis', s=50, alpha=0.7, 
                                      vmin=0, vmax=threshold, label=f'Normal spots ({len(normal_indices)})')
        
        # Draw outlier spots
        if len(outlier_indices) > 0:
            outlier_spots = spots[outlier_indices]
            outlier_distances = distances[outlier_indices]
            scatter_outlier = ax.scatter(outlier_spots[:, 0], outlier_spots[:, 1], 
                                       c='red', s=80, alpha=0.8, marker='^', 
                                       label=f'Outlier spots ({len(outlier_indices)})')
        
        # Draw skeleton
        if len(skeleton) > 0:
            ax.scatter(skeleton[:, 0], skeleton[:, 1], c='gray', s=10, alpha=0.5, label='Skeleton')
        
        # Add colorbar
        if len(normal_indices) > 0:
            cbar = fig.colorbar(scatter_normal, ax=ax, label='Distance to skeleton (μm)')
            cbar.ax.axhline(y=threshold, color='red', linestyle='--', linewidth=2)
        
        ax.set_xlabel('X (μm)')
        ax.set_ylabel('Y (μm)')
        
        flip_status = ""
        if self.use_y_flip.get():
            flip_status += " (Y-flipped)"
        if self.use_z_flip.get():
            flip_status += " (Z-flipped)"
        
        ax.set_title(f'{self.channel.upper()} - {image_name} - {cell_name} - Outlier XY Projection{flip_status}\nThreshold: {threshold} μm')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axis('equal')
        
        canvas = FigureCanvasTkAgg(fig, parent_frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas.draw()

    def create_outlier_distance_distribution_in_window(self, parent_frame, image_name, cell_name, cell_data, threshold):
        """创建outlier分析的距离分布图"""
        fig = Figure(figsize=(10, 8), dpi=100)
        ax = fig.add_subplot(111)
        
        distances = cell_data['distances']
        outlier_indices = cell_data['outlier_indices']
        
        if len(distances) > 0:
            # 绘制直方图
            ax.hist(distances, bins=20, alpha=0.7, edgecolor='black', color='skyblue')
            
            # 添加统计线
            ax.axvline(np.mean(distances), color='blue', linestyle='--', 
                      label=f'Mean: {np.mean(distances):.3f} μm')
            ax.axvline(np.median(distances), color='green', linestyle='--', 
                      label=f'Median: {np.median(distances):.3f} μm')
            ax.axvline(threshold, color='red', linestyle='-', linewidth=2, 
                      label=f'Threshold: {threshold} μm')
            
            # 高亮outlier区域
            ax.axvspan(threshold, distances.max(), alpha=0.2, color='red', label='Outlier region')
        
        ax.set_xlabel('Distance to skeleton (μm)')
        ax.set_ylabel('Number of molecules')
        
        flip_status = ""
        if self.use_y_flip.get():
            flip_status += " (Y-flipped)"
        if self.use_z_flip.get():
            flip_status += " (Z-flipped)"
        
        ax.set_title(f'{self.channel.upper()} - {image_name} - {cell_name} - Outlier Distance Distribution{flip_status}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        if len(distances) > 0:
            outlier_count = len(outlier_indices)
            total_spots = len(distances)
            normal_count = total_spots - outlier_count
            outlier_ratio = outlier_count / total_spots * 100 if total_spots > 0 else 0
            
            stats_text = f"""Statistics:
Total spots: {total_spots}
Normal spots: {normal_count}
Outlier spots: {outlier_count} ({outlier_ratio:.1f}%)

Distance Statistics:
Mean: {np.mean(distances):.3f} μm
Median: {np.median(distances):.3f} μm
Std: {np.std(distances):.3f} μm
Min: {np.min(distances):.3f} μm
Max: {np.max(distances):.3f} μm

Threshold: {threshold} μm
Transform: {flip_status if flip_status else 'Original'}"""
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)
        
        canvas = FigureCanvasTkAgg(fig, parent_frame)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        canvas.draw()

def main():
    root = tk.Tk()
    app = Interactive3DBatchVisualizer(root)
    root.mainloop()

if __name__ == "__main__":
    main() 