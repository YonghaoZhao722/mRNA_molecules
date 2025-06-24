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

class Interactive3DBatchVisualizer:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive 3D Cell Batch Visualization")
        self.root.geometry("1500x950")

        # 路径设置
        self.base_dir = r'F:\atp\Y333 ATP6 TIM50'
        self.skeleton_root = os.path.join(self.base_dir, 'extracted_cells')
        # self.spots_root = self.base_dir
        self.channel = 'tim50'
        self.spots_root = os.path.join(self.base_dir, f'{self.channel}_spots')

        # 数据结构
        self.available_images = {}  # {image_name: {cell_name: {spots_data, skeleton_file}}}
        self.selected_image = tk.StringVar()
        self.selected_cell = tk.StringVar()
        self.current_cell = None
        self.current_image = None

        # 其它变量
        self.original_spots = None
        self.original_skeleton = None
        self.current_skeleton = None
        self.current_spots = None
        self.distances = None
        self.current_skeleton_translation = np.array([0, 0, 0])  # 当前skeleton的平移量
        self.rotation_x = tk.DoubleVar(value=0.0)
        self.rotation_y = tk.DoubleVar(value=0.0)
        self.rotation_z = tk.DoubleVar(value=0.0)
        self.use_y_flip = tk.BooleanVar(value=True)
        self.use_z_flip = tk.BooleanVar(value=False)
        self.auto_compare_yflip = tk.BooleanVar(value=True)
        self.auto_translate_skeleton = tk.BooleanVar(value=False)

        # 创建界面
        self.create_interface()
        # 扫描所有图像和细胞
        self.scan_available_images_and_cells()
        # 自动加载第一个图像和cell
        if self.available_images:
            first_image = list(self.available_images.keys())[0]
            self.selected_image.set(first_image)
            self.on_image_changed()

    def scan_available_images_and_cells(self):
        """扫描所有图像和细胞，spot文件与图像文件夹通过index和field of view双重匹配"""
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
            # 解析index和fov
            parts = image_folder.split('_')
            # 假设最后两个分别为index和fov
            if len(parts) < 2:
                continue
            index = parts[-2]
            fov = parts[-1]
            # 匹配spot文件：包含 _index_ 且 _fov_，且有_spots
            matched_spot = None
            for spot_file in spot_files:
                if f'_{index}_' in spot_file and f'_{fov}_' in spot_file and '_spots' in spot_file:
                    matched_spot = os.path.join(self.spots_root, spot_file)
                    break
            if not matched_spot:
                print(f"Spot file not found for {image_folder} (index={index}, fov={fov})")
                continue
            # 解析spot文件
            cells_data = self.parse_spots_file(matched_spot)
            # 扫描skeleton文件
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
                        image_cells[corresponding_cell_name] = {
                            'spots_data': cells_data[corresponding_cell_name],
                            'skeleton_file': os.path.join(image_path, skel_file)
                        }
                except Exception as e:
                    print(f"Failed to parse skeleton filename {skel_file}: {e}")
            if image_cells:
                self.available_images[image_folder] = image_cells
        # 更新图像下拉框
        if hasattr(self, 'image_combo'):
            self.image_combo['values'] = list(self.available_images.keys())
            if self.available_images:
                self.selected_image.set(list(self.available_images.keys())[0])

    def parse_spots_file(self, spots_file):
        """与原版一致"""
        with open(spots_file, 'r') as f:
            lines = f.readlines()
        cells_data = {}
        current_cell = None
        reading_spots = False
        for i, line in enumerate(lines):
            line = line.strip()
            if line.startswith('CELL'):
                current_cell = line.split('\t')[1]
                cells_data[current_cell] = {'spots': []}
                reading_spots = False
            elif line == 'SPOTS':
                reading_spots = True
                continue
            elif reading_spots and line and not line.startswith('Pos_Y'):
                if line.startswith('CELL'):
                    current_cell = line.split('\t')[1]
                    cells_data[current_cell] = {'spots': []}
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

    def read_skeleton_data(self, skeleton_file):
        df = pd.read_csv(skeleton_file, sep='\t')
        coords = df[['x', 'y', 'z']].values
        return coords

    def convert_and_scale_coordinates(self, spots_nm, skeleton_coords, voxel_size=(300, 160, 160)):
        """使用正确的坐标转换和缩放逻辑，现在支持Y轴和Z轴翻转，并增加数据格式检查"""
        # 检查spots数据格式
        if spots_nm.ndim != 2 or spots_nm.shape[1] != 3:
            raise ValueError(f"Invalid spots data format: expected 2D array with 3 columns, got shape {spots_nm.shape}")
        
        # 直接从纳米转换为微米（修正逻辑）
        spots_um = spots_nm / 1000.0
        
        # 重新排列为(x, y, z)格式以匹配skeleton
        spots_um_xyz = spots_um[:, [2, 1, 0]]  # 从(z, y, x)到(x, y, z)
        
        # 计算缩放因子（恢复正确的缩放逻辑）
        if len(spots_um_xyz) > 1:
            spots_x_range = spots_um_xyz[:, 0].max() - spots_um_xyz[:, 0].min()
            spots_y_range = spots_um_xyz[:, 1].max() - spots_um_xyz[:, 1].min()
            skeleton_x_range = skeleton_coords[:, 0].max() - skeleton_coords[:, 0].min()
            skeleton_y_range = skeleton_coords[:, 1].max() - skeleton_coords[:, 1].min()
            if spots_x_range > 0 and spots_y_range > 0:
                scale_x = skeleton_x_range / spots_x_range
                scale_y = skeleton_y_range / spots_y_range
            else:
                scale_x = scale_y = 1.0
        else:
            scale_x = scale_y = 1.0
        
        # 创建所有可能的变换组合
        transform_results = {}
        
        # 基础坐标（无翻转）
        base_coords = spots_um_xyz.copy()
        
        # 应用Y轴翻转（如果需要）
        if self.use_y_flip.get():
            base_coords[:, 1] = -base_coords[:, 1]
        
        # 应用Z轴翻转（如果需要）- 围绕细胞范围中心翻转
        if self.use_z_flip.get():
            # 计算当前spots z坐标范围中心
            z_center = (base_coords[:, 2].max() + base_coords[:, 2].min()) / 2
            # 围绕中心翻转
            base_coords[:, 2] = 2 * z_center - base_coords[:, 2]
        
        # 应用缩放和对齐
        spots_scaled = base_coords.copy()
        if len(spots_um_xyz) > 1:
            spots_scaled[:, 0] = (spots_scaled[:, 0] - spots_scaled[:, 0].min()) * scale_x + skeleton_coords[:, 0].min()
            spots_scaled[:, 1] = (spots_scaled[:, 1] - spots_scaled[:, 1].min()) * scale_y + skeleton_coords[:, 1].min()
        else:
            spots_scaled[:, 0] = (skeleton_coords[:, 0].min() + skeleton_coords[:, 0].max()) / 2
            spots_scaled[:, 1] = (skeleton_coords[:, 1].min() + skeleton_coords[:, 1].max()) / 2
        spots_scaled[:, 2] = base_coords[:, 2]
        
        # 为了兼容性，仍然返回original和Y flip版本
        # 但现在这些都受当前设置影响
        return {
            'original': spots_scaled,  # 实际上现在包括所有当前变换
            'flip_y': spots_scaled,    # 保持兼容性
            'scale_x': scale_x,
            'scale_y': scale_y
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
        
        # 添加分隔线
        ttk.Separator(transform_frame, orient='horizontal').pack(fill=tk.X, pady=(10, 10))
        
        # 批量分析选项
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
            # 重置平移信息，因为坐标变换发生了变化
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
        """加载指定细胞数据，增加数据格式检查，即使没有spots也显示skeleton"""
        if image_name not in self.available_images or cell_name not in self.available_images[image_name]:
            print(f"Cell {cell_name} in image {image_name} not available")
            return
        
        try:
            self.current_image = image_name
            self.current_cell = cell_name
            cell_info = self.available_images[image_name][cell_name]
            
            # 加载skeleton数据（总是加载）
            skeleton_coords = self.read_skeleton_data(cell_info['skeleton_file'])
            self.original_skeleton = skeleton_coords.copy()
            self.current_skeleton = skeleton_coords.copy()
            
            # 获取spots坐标
            spots_nm = np.array(cell_info['spots_data']['spots'])
            
            # 检查spots数据格式
            if spots_nm.ndim != 2 or spots_nm.shape[1] != 3 or spots_nm.shape[0] == 0:
                print(f"Warning: Invalid or empty spots data for {image_name}-{cell_name}: shape {spots_nm.shape}")
                # 没有有效spots数据，但仍然显示skeleton
                self.current_spots = None
                self.original_spots = None
            else:
                # 使用正确的坐标转换
                conversion_results = self.convert_and_scale_coordinates(spots_nm, skeleton_coords)
                
                # 现在坐标转换函数处理所有变换选项
                self.current_spots = conversion_results['original'].copy()
                self.original_spots = conversion_results['original'].copy()
            
            # 更新窗口标题
            self.root.title(f"Interactive 3D Cell Batch Visualization - {image_name} - {cell_name}")
            
            spots_count = len(self.original_spots) if self.original_spots is not None else 0
            print(f"Successfully loaded {image_name}-{cell_name} data:")
            print(f"  Spots count: {spots_count}")
            print(f"  Skeleton points count: {len(self.original_skeleton)}")
            
            # 重置旋转
            self.reset_rotation()
            
        except Exception as e:
            print(f"Failed to load {image_name}-{cell_name} data: {e}")
            import traceback
            traceback.print_exc()

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
        自动寻找最优的skeleton平移位置，使得到spots的距离中位数最小
        
        Args:
            spots: spots坐标 (N, 3)
            skeleton: skeleton坐标 (M, 3)  
            search_range: 在每个方向上的搜索范围 (μm)
            step_size: 搜索步长 (μm)
        
        Returns:
            dict: {'optimal_skeleton': 最优平移后的skeleton, 'translation': 最优平移向量, 'median_distance': 最小中位数距离}
        """
        if len(spots) == 0 or len(skeleton) == 0:
            return {'optimal_skeleton': skeleton, 'translation': np.array([0, 0, 0]), 'median_distance': float('inf')}
        
        from scipy.optimize import minimize
        
        def objective_function(translation):
            """目标函数：返回平移后skeleton到spots的距离中位数"""
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
            print(f"优化失败，使用原始skeleton: {e}")
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
        
        # 如果启用了自动平移，先进行平移优化
        if self.auto_translate_skeleton.get():
            translation_result = self.optimize_skeleton_translation(self.current_spots, self.current_skeleton)
            self.current_skeleton = translation_result['optimal_skeleton']
            self.current_skeleton_translation = translation_result['translation']
            print(f"自动平移优化完成，平移量: X={self.current_skeleton_translation[0]:.3f}, Y={self.current_skeleton_translation[1]:.3f}, Z={self.current_skeleton_translation[2]:.3f}")
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
        if len(self.current_spots) > 0:
            scatter = ax.scatter(self.current_spots[:, 0], self.current_spots[:, 1], c=self.distances, cmap='viridis', s=50, alpha=0.7)
            fig.colorbar(scatter, ax=ax, label='Distance to skeleton (μm)')
        if len(self.current_skeleton) > 0:
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
        if len(self.current_skeleton) > 0:
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
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
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
        """批量分析所有图像所有细胞的距离分布，输出聚合直方图和boxplot，并统计阈值占比，带进度条"""
        threshold = simpledialog.askfloat("阈值输入", "请输入距离阈值 (μm)：", minvalue=0.0, initialvalue=0.6)
        if threshold is None:
            return
        # 创建进度条窗口
        progress_win = tk.Toplevel(self.root)
        progress_win.title("批量分析进度")
        progress_label = ttk.Label(progress_win, text="正在分析...")
        progress_label.pack(padx=20, pady=10)
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_win, variable=progress_var, maximum=1.0, length=300)
        progress_bar.pack(padx=20, pady=10)
        self.root.update()
        all_distances = []
        summary_rows = []
        cell_spot_counts = []
        total = sum(len(cells) for cells in self.available_images.values())
        done = 0
        for image_name, cells in self.available_images.items():
            for cell_name, cell_info in cells.items():
                spots_nm = np.array(cell_info['spots_data']['spots'])
                skeleton_coords = self.read_skeleton_data(cell_info['skeleton_file'])
                # 跳过无spot或无skeleton的cell
                if spots_nm.ndim != 2 or spots_nm.shape[0] == 0 or skeleton_coords.shape[0] == 0:
                    done += 1
                    progress_var.set(done / total)
                    self.root.update()
                    continue
                try:
                    # 根据设置决定是否对比Y-flip
                    if self.auto_compare_yflip.get():
                        # 自动对比模式：分别计算y_flip True/False，选择最优
                        best_result = None
                        best_mean = None
                        for y_flip in [True, False]:
                            self.use_y_flip.set(y_flip)
                            conversion_results = self.convert_and_scale_coordinates(spots_nm, skeleton_coords)
                            spots = conversion_results['original']
                            skeleton = skeleton_coords
                            
                            # 应用自动平移优化（如果启用）
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
                        # 用最优结果
                        spots = best_result['spots']
                        skeleton = best_result['skeleton']
                        min_distances = best_result['min_distances']
                        y_flip_used = best_result['y_flip']
                        skeleton_translation = best_result['skeleton_translation']
                    else:
                        # 使用当前界面设置，不进行对比
                        conversion_results = self.convert_and_scale_coordinates(spots_nm, skeleton_coords)
                        spots = conversion_results['original']
                        skeleton = skeleton_coords
                        distances = cdist(spots, skeleton)
                        min_distances = np.min(distances, axis=1)
                        y_flip_used = self.use_y_flip.get()
                    
                    # 应用自动平移优化（如果启用）
                    skeleton_translation = np.array([0, 0, 0])  # 默认无平移
                    if self.auto_translate_skeleton.get():
                        translation_result = self.optimize_skeleton_translation(spots, skeleton)
                        skeleton = translation_result['optimal_skeleton']
                        skeleton_translation = translation_result['translation']
                        # 重新计算距离
                        distances = cdist(spots, skeleton)
                        min_distances = np.min(distances, axis=1)
                    below = (min_distances < threshold)
                    for d in min_distances:
                        all_distances.append(d)
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
        # 保存所有距离分布
        output_dir = os.path.join(self.base_dir, 'interactive_batch_results')
        os.makedirs(output_dir, exist_ok=True)
        import pandas as pd
        dist_df = pd.DataFrame({'distance': all_distances})
        dist_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_distances.csv'), index=False)
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(os.path.join(output_dir, f'{self.channel}_all_cells_summary.csv'), index=False)
        # 画聚合直方图和boxplot
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
            # 画mean、median、threshold竖线
            plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
            plt.axvline(median_val, color='orange', linestyle='--', label=f'Median: {median_val:.2f}')
            plt.axvline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            # 标注低于阈值的占比、细胞总数、平均spot数
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
            # 画threshold线
            plt.axhline(threshold, color='purple', linestyle='-', label=f'Threshold: {threshold:.2f}')
            # 标注低于阈值的占比、细胞总数、平均spot数
            plt.text(1.05, threshold, f'Ratio < threshold: {below_ratio*100:.1f}%\nNumber of cells: {cell_count}\nMean spots per cell: {mean_spots_per_cell:.1f}',
                     transform=plt.gca().get_yaxis_transform(which='grid'), ha='left', va='center', fontsize=12,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{self.channel}_all_cells_boxplot.png'))
            plt.close()
        progress_win.destroy()
        messagebox.showinfo("批量分析完成", f"所有细胞距离分布直方图和boxplot已保存到：\n{output_dir}")

def main():
    root = tk.Tk()
    app = Interactive3DBatchVisualizer(root)
    root.mainloop()

if __name__ == "__main__":
    main() 