#!/usr/bin/env python3
"""
启动Cell 2交互式3D可视化工具

使用方法:
    python run_interactive_cell2.py

功能:
- 加载Cell 2的spots和skeleton数据
- 提供3D交互式可视化
- 允许用户调整skeleton的方向（X/Y/Z轴旋转）
- 计算并显示XY投影和距离分布
- 保存分析结果
"""

import sys
import os

# 添加当前目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from interactive_3d_cell2_visualization import main
    
    if __name__ == "__main__":
        print("启动Cell 2交互式3D可视化工具...")
        print("=" * 50)
        print("使用说明:")
        print("1. 程序会自动加载Cell 2的数据")
        print("2. 使用左侧滑条调整skeleton的方向")
        print("3. 点击'计算距离分析'按钮查看结果")
        print("4. 点击'保存结果'按钮保存分析数据")
        print("=" * 50)
        
        main()
        
except ImportError as e:
    print(f"导入错误: {e}")
    print("请确保所有依赖包已安装:")
    print("pip install -r requirements.txt")
except Exception as e:
    print(f"启动失败: {e}")
    print("请检查数据文件是否存在和路径是否正确") 