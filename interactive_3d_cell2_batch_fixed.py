import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import threading
import tkinter.simpledialog as simpledialog
import json

class Interactive3DBatchVisualizerFixed:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive 3D Cell Batch Visualization (Fixed)")
        self.root.geometry("1500x950")

        # Path settings
        self.base_dir = r'Y333 ATP6 ATP2'
        self.skeleton_root = os.path.join(self.base_dir, 'extracted_cells')
        self.channel = 'atp6'
        self.spots_root = os.path.join(self.base_dir, f'{self.channel}_spots')

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
        self.distances = None
        self.current_skeleton_translation = np.array([0, 0, 0])
        
        # Control variables
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

    def load_spots_fishquant_analyze_method(self, file_path, cell_number=1, flip_y=True, mapping_data=None, silent=False):
        """
        完全复制analyze_alignment.py中load_spots_fishquant的正确方法
        修复：正确处理某些细胞没有spots数据的情况
        """
        # 读取文件找到SPOTS部分
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # 找到像素大小信息
        pixel_xy, pixel_z = None, None
        for i, line in enumerate(lines):
            if line.startswith('Pix-XY'):
                # 像素大小数值在下一行
                if i + 1 < len(lines):
                    next_line = lines[i + 1].strip()
                    parts = next_line.split('\t')
                    if len(parts) >= 2:
                        pixel_xy = float(parts[0])
                        pixel_z = float(parts[1])
                break
        
        if not silent:
            print(f"FISH-QUANT像素大小: XY={pixel_xy}nm, Z={pixel_z}nm")
        
        # 找到指定细胞的SPOTS数据
        target_cell = f"CELL\tCell_{cell_number}"
        cell_found = False
        spots_start = -1
        spots_end = -1
        
        for i, line in enumerate(lines):
            line_stripped = line.strip()
            
            # 遇到目标细胞标记
            if line_stripped == target_cell:
                cell_found = True
                if not silent:
                    print(f"找到目标细胞: Cell_{cell_number}")
                continue
            
            # 如果已找到目标细胞，且遇到SPOTS标记（说明这个细胞有spots数据）
            if cell_found and line.startswith('Pos_Y'):
                spots_start = i
                continue
            
            # 如果已找到目标细胞和spots开始位置，遇到下一个CELL标记表示结束
            if cell_found and spots_start != -1 and line.startswith('CELL'):
                spots_end = i
                break
                
            # 关键修复：如果已找到目标细胞但还没找到spots开始位置，
            # 遇到下一个CELL标记，说明目标细胞没有spots数据
            if cell_found and spots_start == -1 and line.startswith('CELL'):
                if not silent:
                    print(f"Cell_{cell_number} 没有spots数据（在找到spots开始位置前遇到了下一个细胞）")
                raise ValueError(f"Cell_{cell_number} 没有spots数据")
        
        # 如果没找到下一个CELL标记，说明是最后一个细胞
        if cell_found and spots_start != -1 and spots_end == -1:
            spots_end = len(lines)
        
        if not cell_found:
            raise ValueError(f"未找到Cell_{cell_number}的数据")
        
        if spots_start == -1:
            raise ValueError(f"Cell_{cell_number} 没有spots数据")
        
        # 读取该细胞的spots数据
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
        
        # 如果需要翻转Y轴
        if flip_y and mapping_data:
            # 使用mapping数据中的crop区域信息进行精确翻转
            cell_name = f"cell_{cell_number:03d}"
            if cell_name in mapping_data:
                crop_info = mapping_data[cell_name]['crop_region']
                y_start = crop_info['y_start']  # 像素坐标
                y_end = crop_info['y_end']      # 像素坐标
                
                # 转换为纳米坐标 (像素 * 像素大小)
                y_start_nm = y_start * pixel_xy
                y_end_nm = y_end * pixel_xy
                
                # 基于crop区域进行翻转：new_y = (y_start + y_end) - old_y
                flip_center = y_start_nm + y_end_nm
                coords[:, 0] = flip_center - coords[:, 0]
                
                if not silent:
                    print(f"Y轴翻转: 基于crop区域 y_start={y_start_nm:.1f}nm, y_end={y_end_nm:.1f}nm")
            else:
                if not silent:
                    print(f"警告: 未找到{cell_name}的mapping信息，跳过Y轴翻转")
        elif flip_y:
            if not silent:
                print("警告: 需要mapping数据进行精确Y轴翻转，使用简单翻转")
            # 简单翻转：基于当前细胞spots的范围
            if len(coords) > 0:
                y_min, y_max = coords[:, 0].min(), coords[:, 0].max()
                flip_center = y_min + y_max
                coords[:, 0] = flip_center - coords[:, 0]
                if not silent:
                    print(f"Y轴翻转: 基于细胞范围，翻转中心={flip_center/2:.1f}nm")
        
        if not silent:
            print(f"Cell_{cell_number} Spots数据加载完成: {len(coords)} 个点")
            print(f"坐标范围: X=[{coords[:,1].min():.1f}, {coords[:,1].max():.1f}], Y=[{coords[:,0].min():.1f}, {coords[:,0].max():.1f}], Z=[{coords[:,2].min():.1f}, {coords[:,2].max():.1f}]")
        
        return coords, pixel_xy, pixel_z

    def load_skeleton_txt_analyze_method(self, file_path, mapping_data=None, cell_name="cell_001", pixel_size_xy=0.0645, silent=False):
        """
        完全复制analyze_alignment.py中load_skeleton_txt的正确方法
        """
        df = pd.read_csv(file_path, sep='\t')
        # 提取坐标列 (x, y, z) - 这些是相对于截取区域的坐标
        coords = df[['x', 'y', 'z']].values
        if not silent:
            print(f"Skeleton数据加载完成: {len(coords)} 个点")
            print(f"相对坐标范围: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
        
        # 如果提供了mapping数据，转换为绝对坐标
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
                
                if not silent:
                    print(f"应用坐标偏移: X+{offset_x_um:.3f}μm (像素{x_offset}), Y+{offset_y_um:.3f}μm (像素{y_offset})")
                    print(f"绝对坐标范围: X=[{coords[:,0].min():.3f}, {coords[:,0].max():.3f}], Y=[{coords[:,1].min():.3f}, {coords[:,1].max():.3f}], Z=[{coords[:,2].min():.3f}, {coords[:,2].max():.3f}]")
            else:
                if not silent:
                    print(f"警告: 在mapping文件中未找到{cell_name}")
        else:
            if not silent:
                if mapping_data:
                    print(f"警告: mapping文件不存在")
                print("使用相对坐标（未应用偏移）")
        
        return coords 

    def scan_available_images_and_cells(self):
        """扫描所有图像和细胞，通过索引和视野双重匹配spot文件与图像文件夹
        修复：准确检测哪些细胞真正有spots数据"""
        self.available_images = {}
        if not os.path.exists(self.skeleton_root):
            print(f"Skeleton root not found: {self.skeleton_root}")
            return
            
        # 获取所有spot文件
        spot_files = [f for f in os.listdir(self.spots_root) if f.endswith('.txt') and '_spots' in f]
        
        for image_folder in os.listdir(self.skeleton_root):
            image_path = os.path.join(self.skeleton_root, image_folder)
            if not os.path.isdir(image_path):
                continue

            # 加载坐标映射
            mapping_file = os.path.join(image_path, 'coordinate_mapping.json')
            if os.path.exists(mapping_file):
                with open(mapping_file, 'r') as f:
                    self.coordinate_mappings[image_folder] = json.load(f)
            else:
                print(f"Warning: coordinate_mapping.json not found in {image_path}")
                continue

            # 解析索引和视野
            parts = image_folder.split('_')
            if len(parts) < 2:
                continue
            index = parts[-2]
            fov = parts[-1]
            
            # 匹配spot文件
            matched_spot = None
            for spot_file in spot_files:
                if f'_{index}_' in spot_file and f'_{fov}_' in spot_file and '_spots' in spot_file:
                    matched_spot = os.path.join(self.spots_root, spot_file)
                    break
            
            if not matched_spot:
                print(f"Spot file not found for {image_folder} (index={index}, fov={fov})")
                continue

            # 扫描skeleton文件
            skeleton_files = [f for f in os.listdir(image_path) if f.startswith('cell_') and f.endswith('.txt')]
            image_cells = {}
            
            for skel_file in skeleton_files:
                try:
                    cell_num = int(skel_file.split('_')[1].split('.')[0])
                    cell_name = f"Cell_{cell_num}"
                    
                    # 严格检查这个细胞是否在spot文件中有实际的spots数据
                    try:
                        # 使用修复后的函数来严格验证spots数据
                        test_coords, _, _ = self.load_spots_fishquant_analyze_method(
                            matched_spot, 
                            cell_number=cell_num, 
                            flip_y=False,  # 测试时不翻转
                            mapping_data=None,  # 测试时不用mapping
                            silent=True  # 静默模式
                        )
                        
                        # 如果能成功加载且有数据，则添加到可用列表
                        if len(test_coords) > 0:
                            image_cells[cell_name] = {
                                'spots_file': matched_spot,
                                'skeleton_file': os.path.join(image_path, skel_file)
                            }
                            print(f"Found {image_folder}-{cell_name} (verified with {len(test_coords)} spots)")
                        else:
                            print(f"Skipping {image_folder}-{cell_name}: spots data is empty")
                            
                    except ValueError as e:
                        # 如果抛出"没有spots数据"的异常，跳过这个细胞
                        print(f"Skipping {image_folder}-{cell_name}: {str(e)}")
                    except Exception as e:
                        print(f"Error checking spots file for {cell_name}: {e}")
                        
                except Exception as e:
                    print(f"Failed to parse skeleton filename {skel_file}: {e}")
            
            if image_cells:
                self.available_images[image_folder] = image_cells
        
        # 更新图像下拉框
        if hasattr(self, 'image_combo'):
            self.image_combo['values'] = list(self.available_images.keys())
            if self.available_images:
                self.selected_image.set(list(self.available_images.keys())[0])

    def load_cell_data_analyze_method(self, image_name, cell_name):
        """
        使用analyze_alignment.py验证过的正确方法加载细胞数据
        """
        if image_name not in self.available_images or cell_name not in self.available_images[image_name]:
            print(f"Cell {cell_name} in image {image_name} not available")
            return
        
        try:
            self.current_image = image_name
            self.current_cell = cell_name
            cell_info = self.available_images[image_name][cell_name]
            
            # 获取细胞编号和映射数据
            cell_number = int(cell_name.split("_")[1])
            cell_id_str = f"{cell_number:03d}"
            mapping_cell_name = f'cell_{cell_id_str}'
            mapping_data = self.coordinate_mappings.get(image_name, {})
            
            print(f"=== 加载 {image_name} - {cell_name} ===")
            
            # 1. 加载spots数据（使用analyze_alignment.py的方法）
            spots_coords, pixel_xy, pixel_z = self.load_spots_fishquant_analyze_method(
                cell_info['spots_file'], 
                cell_number=cell_number, 
                flip_y=self.use_y_flip.get(), 
                mapping_data=mapping_data
            )
            
            # 2. 加载skeleton数据（使用analyze_alignment.py的方法）
            skeleton_coords = self.load_skeleton_txt_analyze_method(
                cell_info['skeleton_file'], 
                mapping_data=mapping_data,
                cell_name=mapping_cell_name,
                pixel_size_xy=pixel_xy/1000  # 转换为微米
            )
            
            # 3. 坐标转换（完全按照analyze_alignment.py的逻辑）
            # 重新排列为(x, y, z)格式以匹配skeleton
            spots_nm_xyz = spots_coords[:, [1, 0, 2]]  # from (y, x, z) to (x, y, z)
            
            # 转换为微米单位
            spots_um_xyz = spots_nm_xyz / 1000.0
            
            print(f"转换后spots坐标范围: X=[{spots_um_xyz[:,0].min():.3f}, {spots_um_xyz[:,0].max():.3f}], Y=[{spots_um_xyz[:,1].min():.3f}, {spots_um_xyz[:,1].max():.3f}], Z=[{spots_um_xyz[:,2].min():.3f}, {spots_um_xyz[:,2].max():.3f}]")
            
            # 4. 中心校正分析（完全复制analyze_alignment.py的逻辑）
            skeleton_center = np.mean(skeleton_coords, axis=0)
            spots_center = np.mean(spots_um_xyz, axis=0)
            center_diff = skeleton_center - spots_center
            center_distance = np.linalg.norm(center_diff)
            
            print(f"中心校正分析:")
            print(f"  Skeleton中心: X={skeleton_center[0]:.3f}, Y={skeleton_center[1]:.3f}, Z={skeleton_center[2]:.3f}")
            print(f"  Spots中心: X={spots_center[0]:.3f}, Y={spots_center[1]:.3f}, Z={spots_center[2]:.3f}")
            print(f"  中心距离: {center_distance:.3f}μm")
            
            # 应用中心校正
            spots_corrected = spots_um_xyz - center_diff
            
            # 保存两个版本的数据以匹配analyze_alignment.py的行为
            self.spots_before_correction = spots_um_xyz.copy()  # 中心校正前（用于图表显示，匹配analyze_alignment.py）
            self.spots_after_correction = spots_corrected.copy()  # 中心校正后（用于最终统计）
            
            # 默认使用中心校正前的数据进行可视化（与analyze_alignment.py的图表一致）
            self.current_spots = self.spots_before_correction.copy()
            self.original_spots = self.spots_before_correction.copy()
            self.current_skeleton = skeleton_coords.copy()
            self.original_skeleton = skeleton_coords.copy()
            
            print(f"校正前spots坐标范围（用于可视化）: X=[{self.current_spots[:,0].min():.3f}, {self.current_spots[:,0].max():.3f}], Y=[{self.current_spots[:,1].min():.3f}, {self.current_spots[:,1].max():.3f}], Z=[{self.current_spots[:,2].min():.3f}, {self.current_spots[:,2].max():.3f}]")
            print(f"校正后spots坐标范围: X=[{spots_corrected[:,0].min():.3f}, {spots_corrected[:,0].max():.3f}], Y=[{spots_corrected[:,1].min():.3f}, {spots_corrected[:,1].max():.3f}], Z=[{spots_corrected[:,2].min():.3f}, {spots_corrected[:,2].max():.3f}]")
            
            # 更新窗口标题
            self.root.title(f"Interactive 3D Cell Batch Visualization (Fixed) - {image_name} - {cell_name}")
            
            print(f"Successfully loaded {image_name}-{cell_name} data:")
            print(f"  Spots count: {len(self.original_spots)}")
            print(f"  Skeleton points count: {len(self.original_skeleton)}")
            
            # 重置旋转
            self.reset_rotation()
            
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
        self.create_rotation_controls(control_frame)
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
        
        z_flip_check = ttk.Checkbutton(transform_frame, text="Z-axis Flip", variable=self.use_z_flip, command=self.on_transform_change)
        z_flip_check.pack(anchor=tk.W)
        
        # Add separator
        ttk.Separator(transform_frame, orient='horizontal').pack(fill=tk.X, pady=(10, 10))
        
        # Batch analysis options
        ttk.Label(transform_frame, text="Batch Analysis Options:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        auto_compare_check = ttk.Checkbutton(transform_frame, text="Auto Y-flip comparison", variable=self.auto_compare_yflip)
        auto_compare_check.pack(anchor=tk.W)
        auto_translate_check = ttk.Checkbutton(transform_frame, text="Use skeleton auto-translation", variable=self.auto_translate_skeleton)
        auto_translate_check.pack(anchor=tk.W)

    def create_rotation_controls(self, parent):
        rotation_frame = ttk.LabelFrame(parent, text="Skeleton Rotation", padding=10)
        rotation_frame.pack(fill=tk.X, pady=(0, 15))
        
        # X轴旋转
        x_frame = ttk.Frame(rotation_frame)
        x_frame.pack(fill=tk.X, pady=2)
        ttk.Label(x_frame, text="X-axis:", width=6).pack(side=tk.LEFT)
        x_scale = ttk.Scale(x_frame, from_=-180, to=180, variable=self.rotation_x, orient=tk.HORIZONTAL, command=self.on_rotation_change)
        x_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        self.x_value_label = ttk.Label(x_frame, text="0°", width=6)
        self.x_value_label.pack(side=tk.RIGHT)
        
        # Y轴旋转
        y_frame = ttk.Frame(rotation_frame)
        y_frame.pack(fill=tk.X, pady=2)
        ttk.Label(y_frame, text="Y-axis:", width=6).pack(side=tk.LEFT)
        y_scale = ttk.Scale(y_frame, from_=-180, to=180, variable=self.rotation_y, orient=tk.HORIZONTAL, command=self.on_rotation_change)
        y_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        self.y_value_label = ttk.Label(y_frame, text="0°", width=6)
        self.y_value_label.pack(side=tk.RIGHT)
        
        # Z轴旋转
        z_frame = ttk.Frame(rotation_frame)
        z_frame.pack(fill=tk.X, pady=2)
        ttk.Label(z_frame, text="Z-axis:", width=6).pack(side=tk.LEFT)
        z_scale = ttk.Scale(z_frame, from_=-180, to=180, variable=self.rotation_z, orient=tk.HORIZONTAL, command=self.on_rotation_change)
        z_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        self.z_value_label = ttk.Label(z_frame, text="0°", width=6)
        self.z_value_label.pack(side=tk.RIGHT)

    def create_action_buttons(self, parent):
        button_frame = ttk.LabelFrame(parent, text="Actions", padding=10)
        button_frame.pack(fill=tk.X, pady=(0, 15))
        
        reset_btn = ttk.Button(button_frame, text="Reset Rotation", command=self.reset_rotation)
        reset_btn.pack(fill=tk.X, pady=2)
        
        calculate_btn = ttk.Button(button_frame, text="Calculate Distance", command=self.calculate_and_show_results)
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
            self.ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], self.current_spots[:, 2], 
                          c='blue', s=30, alpha=0.6, label=f'Spots ({len(self.current_spots)})')
        
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
            self.ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], self.current_skeleton[:, 2], 
                          c='red', s=10, alpha=0.8, label=f'Skeleton ({len(self.current_skeleton)})')
        
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

    def calculate_distances_analyze_method(self):
        """
        完全复制analyze_alignment.py中calculate_alignment_after_correction的距离计算方法
        使用中心校正前的数据和采样，以匹配analyze_alignment.py的图表
        """
        if not hasattr(self, 'spots_before_correction') or self.current_skeleton is None:
            return np.array([])
        
        spots_coords_um = self.spots_before_correction
        skeleton_coords = self.current_skeleton
        
        if len(spots_coords_um) == 0 or len(skeleton_coords) == 0:
            return np.array([])
        
        # 采样逻辑（完全复制analyze_alignment.py）
        if len(skeleton_coords) > 1000:
            skeleton_sample = skeleton_coords[::len(skeleton_coords)//1000]
        else:
            skeleton_sample = skeleton_coords
        
        if len(spots_coords_um) > 500:
            spots_sample = spots_coords_um[::len(spots_coords_um)//500]
        else:
            spots_sample = spots_coords_um
        
        distances = cdist(spots_sample, skeleton_sample)
        min_distances = np.min(distances, axis=1)
        return min_distances

    def calculate_final_distances(self):
        """
        计算最终对齐统计（使用中心校正后的数据）
        """
        if not hasattr(self, 'spots_after_correction') or self.current_skeleton is None:
            return np.array([])
        
        spots_corrected = self.spots_after_correction
        skeleton_coords = self.current_skeleton
        
        if len(spots_corrected) == 0 or len(skeleton_coords) == 0:
            return np.array([])
        
        # 使用和analyze_alignment.py相同的采样策略
        spots_sample = spots_corrected[:500] if len(spots_corrected) > 500 else spots_corrected
        skeleton_sample = skeleton_coords[:1000] if len(skeleton_coords) > 1000 else skeleton_coords
        
        distances = cdist(spots_sample, skeleton_sample)
        min_distances = np.min(distances, axis=1)
        return min_distances

    def calculate_and_show_results(self):
        if self.current_spots is None or self.current_skeleton is None:
            messagebox.showwarning("Warning", "请先选择一个细胞")
            return
        
        # 计算两种距离（匹配analyze_alignment.py的行为）
        self.distances = self.calculate_distances_analyze_method()  # 用于图表显示
        self.final_distances = self.calculate_final_distances()      # 用于最终统计
        
        if len(self.distances) == 0:
            messagebox.showwarning("Warning", "无法计算距离，请检查数据")
            return
            
        self.create_result_window()
        
        print(f"{self.current_image} - {self.current_cell} 计算完成:")
        print(f"  图表显示距离（中心校正前）:")
        print(f"    平均距离: {np.mean(self.distances):.3f} μm")
        print(f"    中位数距离: {np.median(self.distances):.3f} μm")
        print(f"    标准差: {np.std(self.distances):.3f} μm")
        
        if len(self.final_distances) > 0:
            print(f"  最终对齐统计（中心校正后）:")
            print(f"    平均距离: {np.mean(self.final_distances):.3f} μm")
            print(f"    中位数距离: {np.median(self.final_distances):.3f} μm")
            print(f"    标准差: {np.std(self.final_distances):.3f} μm")

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
            ax.scatter(self.current_skeleton[:, 0], self.current_skeleton[:, 1], 
                      c='red', s=10, alpha=0.5, label='Skeleton')
        
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
        
        if self.current_skeleton is not None and len(self.current_skeleton) > 0:
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
            
            # 显示图表数据统计（中心校正前，匹配analyze_alignment.py图表）
            stats_text = f"""Graph Statistics (Pre-Correction):
Count: {len(self.distances)}
Mean: {np.mean(self.distances):.3f} μm
Median: {np.median(self.distances):.3f} μm
Std: {np.std(self.distances):.3f} μm
Min: {np.min(self.distances):.3f} μm
Max: {np.max(self.distances):.3f} μm"""

            # 如果有最终统计数据，也显示
            if hasattr(self, 'final_distances') and len(self.final_distances) > 0:
                stats_text += f"""

Final Statistics (Post-Correction):
Count: {len(self.final_distances)}
Mean: {np.mean(self.final_distances):.3f} μm
Median: {np.median(self.final_distances):.3f} μm
Std: {np.std(self.final_distances):.3f} μm"""

            stats_text += f"""

Transform: {flip_status if flip_status else 'Original'}
Rotation: X={self.rotation_x.get():.0f}° Y={self.rotation_y.get():.0f}° Z={self.rotation_z.get():.0f}°
Note: Graph shows pre-correction data (matches analyze_alignment.py)"""
            
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

    def save_results(self):
        if not hasattr(self, 'distances') or len(self.distances) == 0:
            messagebox.showwarning("Warning", "没有结果可保存，请先计算距离分析")
            return
        if not self.current_cell or not self.current_image:
            messagebox.showwarning("Warning", "没有选择细胞或图像")
            return
        
        try:
            output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
            os.makedirs(output_dir, exist_ok=True)
            
            # 保存详细数据（使用图表显示的中心校正前数据）
            detail_data = pd.DataFrame({
                'spot_x_um': self.spots_before_correction[:, 0],
                'spot_y_um': self.spots_before_correction[:, 1],
                'spot_z_um': self.spots_before_correction[:, 2],
                'distance_to_skeleton_um': self.distances,
                'rotation_x': self.rotation_x.get(),
                'rotation_y': self.rotation_y.get(),
                'rotation_z': self.rotation_z.get(),
                'y_flip_used': self.use_y_flip.get(),
                'z_flip_used': self.use_z_flip.get(),
                'auto_compare_enabled': self.auto_compare_yflip.get(),
                'auto_translate_enabled': self.auto_translate_skeleton.get(),
                'image_name': self.current_image,
                'cell_name': self.current_cell
            })
            detail_file = os.path.join(output_dir, f'{self.current_image}_{self.current_cell}_interactive_analysis_fixed.csv')
            detail_data.to_csv(detail_file, index=False)
            
            # 保存skeleton数据
            skeleton_data = pd.DataFrame({
                'skeleton_x_um': self.current_skeleton[:, 0],
                'skeleton_y_um': self.current_skeleton[:, 1],
                'skeleton_z_um': self.current_skeleton[:, 2],
            })
            skeleton_file = os.path.join(output_dir, f'{self.current_image}_{self.current_cell}_rotated_skeleton_fixed.csv')
            skeleton_data.to_csv(skeleton_file, index=False)
            
            # 保存摘要
            summary_data = {
                'image_name': self.current_image,
                'cell_name': self.current_cell,
                'num_spots': len(self.spots_before_correction),
                'num_skeleton_points': len(self.current_skeleton),
                'y_flip_used': self.use_y_flip.get(),
                'z_flip_used': self.use_z_flip.get(),
                'auto_compare_enabled': self.auto_compare_yflip.get(),
                'auto_translate_enabled': self.auto_translate_skeleton.get(),
                'rotation_x': self.rotation_x.get(),
                'rotation_y': self.rotation_y.get(),
                'rotation_z': self.rotation_z.get(),
                'mean_distance_pre_correction': np.mean(self.distances),
                'median_distance_pre_correction': np.median(self.distances),
                'std_distance_pre_correction': np.std(self.distances),
                'min_distance_pre_correction': np.min(self.distances),
                'max_distance_pre_correction': np.max(self.distances)
            }
            
            # 如果有最终统计数据，也加入
            if hasattr(self, 'final_distances') and len(self.final_distances) > 0:
                summary_data.update({
                    'mean_distance_post_correction': np.mean(self.final_distances),
                    'median_distance_post_correction': np.median(self.final_distances),
                    'std_distance_post_correction': np.std(self.final_distances),
                    'min_distance_post_correction': np.min(self.final_distances),
                    'max_distance_post_correction': np.max(self.final_distances)
                })
            
            summary_df = pd.DataFrame([summary_data])
            summary_file = os.path.join(output_dir, f'{self.current_image}_{self.current_cell}_summary_fixed.csv')
            summary_df.to_csv(summary_file, index=False)
            
            messagebox.showinfo("Success", f"结果已保存到:\n{output_dir}")
        except Exception as e:
            messagebox.showerror("Error", f"保存失败: {str(e)}")

    def batch_distance_analysis(self):
        """批量分析所有图像和细胞的距离分布，使用修复后的正确坐标转换方法"""
        threshold = simpledialog.askfloat("阈值输入", "请输入距离阈值 (μm):", minvalue=0.0, initialvalue=0.6)
        if threshold is None:
            return
        
        # 创建进度条窗口
        progress_win = tk.Toplevel(self.root)
        progress_win.title("批量分析进度")
        progress_label = ttk.Label(progress_win, text="分析中...")
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
                            
                            # 计算距离（使用采样）
                            skeleton_sample = skeleton[:1000] if len(skeleton) > 1000 else skeleton
                            spots_sample = spots[:500] if len(spots) > 500 else spots
                            distances = cdist(spots_sample, skeleton_sample)
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
                        
                        skeleton_sample = skeleton[:1000] if len(skeleton) > 1000 else skeleton
                        spots_sample = spots[:500] if len(spots) > 500 else spots
                        distances = cdist(spots_sample, skeleton_sample)
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
        dist_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_distances_fixed.csv'), index=False)
        
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_summary_fixed.csv'), index=False)
        
        if all_spots_details:
            spots_details_df = pd.DataFrame(all_spots_details)
            spots_details_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_spots_details_fixed.csv'), index=False)
            
            exceeds_threshold_df = spots_details_df[spots_details_df['exceeds_threshold'] == True]
            exceeds_threshold_df.to_csv(os.path.join(output_dir, f'{self.channel}_spots_exceed_threshold_{threshold}_fixed.csv'), index=False)
        
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
            plt.title(f'{self.channel.upper()} - All Cells Distance Distribution (Fixed Method) (N={len(dist_df)})')
            plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
            plt.axvline(median_val, color='orange', linestyle='--', label=f'Median: {median_val:.2f}')
            plt.axvline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            plt.text(0.98, 0.95, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().transAxes, ha='right', va='top', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_histogram_fixed.png'))
            plt.close()
            
            plt.figure(figsize=(8, 6))
            plt.boxplot(dist_df['distance'], vert=True, patch_artist=True, boxprops=dict(facecolor='lightgreen'))
            plt.ylabel('Distance to skeleton (μm)')
            plt.title(f'{self.channel.upper()} - All Cells Distance Boxplot (Fixed Method) (N={len(dist_df)})')
            plt.axhline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            plt.text(1.05, threshold, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().get_yaxis_transform(which='grid'), ha='left', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_boxplot_fixed.png'))
            plt.close()
        
        progress_win.destroy()
        exceed_count = len(exceeds_threshold_df) if 'exceeds_threshold_df' in locals() else 0
        total_spots = len(all_spots_details) if all_spots_details else 0
        messagebox.showinfo("批量分析完成", f"所有细胞的距离分布直方图和箱线图已保存到:\n{output_dir}\n\n总spots: {total_spots}\n超出阈值 {threshold} μm的spots: {exceed_count}")

    def analyse_outliers(self):
        """分析outliers并在新窗口中显示包含outliers的图像和细胞"""
        # 询问阈值
        threshold = simpledialog.askfloat("Outlier阈值输入", "请输入outlier距离阈值 (μm)：", minvalue=0.0, initialvalue=0.6)
        if threshold is None:
            return
        
        # 创建进度条窗口
        progress_win = tk.Toplevel(self.root)
        progress_win.title("Analyse Outliers Progress")
        progress_label = ttk.Label(progress_win, text="分析outliers中...")
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
                try:
                    # 获取细胞编号和映射数据
                    cell_number = int(cell_name.split("_")[1])
                    cell_id_str = f"{cell_number:03d}"
                    mapping_cell_name = f'cell_{cell_id_str}'
                    mapping_data = self.coordinate_mappings.get(image_name, {})
                    
                    # 使用正确的方法加载数据
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
                    
                    # 坐标转换
                    spots_nm_xyz = spots_coords[:, [1, 0, 2]]
                    spots_um_xyz = spots_nm_xyz / 1000.0
                    
                    # 中心校正
                    skeleton_center = np.mean(skeleton_coords, axis=0)
                    spots_center = np.mean(spots_um_xyz, axis=0)
                    center_diff = skeleton_center - spots_center
                    spots_corrected = spots_um_xyz - center_diff
                    
                    spots = spots_um_xyz  # 使用中心校正前的数据进行可视化
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
                    # 静默跳过无法处理的细胞（例如未找到spots数据）
                    pass
                
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
        outlier_window.title(f"Outlier Analysis (Fixed Method) - Threshold: {threshold} μm")
        outlier_window.geometry("1500x950")
        
        # 创建主框架
        main_frame = ttk.Frame(outlier_window)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 控制面板
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # 标题
        title_label = ttk.Label(control_frame, text=f"Outlier Analysis (Fixed)\nThreshold: {threshold} μm", 
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
            
            outlier_ax.set_title(f'{self.channel.upper()} - {image_name} - {cell_name} - Outlier Analysis (Fixed){flip_status}\nThreshold: {threshold_val} μm')
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
        
        # 操作按钮
        action_frame = ttk.LabelFrame(control_frame, text="Actions", padding=10)
        action_frame.pack(fill=tk.X, pady=(15, 0))
        
        def save_outlier_details():
            """保存所有outlier详细信息到CSV"""
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
                outlier_file = os.path.join(output_dir, f'{self.channel}_outlier_analysis_{threshold}_fixed.csv')
                outlier_df.to_csv(outlier_file, index=False)
                messagebox.showinfo("保存成功", f"Outlier详细信息已保存到:\n{outlier_file}")
        
        # 添加按钮
        save_all_btn = ttk.Button(action_frame, text="Save All Outlier Details", command=save_outlier_details)
        save_all_btn.pack(fill=tk.X, pady=2)


def main():
    root = tk.Tk()
    app = Interactive3DBatchVisualizerFixed(root)
    root.mainloop()

if __name__ == "__main__":
    main() 