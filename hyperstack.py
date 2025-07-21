import os
import tifffile as tiff
import numpy as np

def convert_hyperstack_to_zstack(input_path, output_path=None):
    # 读取图像
    img = tiff.imread(input_path)
    img_shape = img.shape
    print(f"Reading {os.path.basename(input_path)} - shape: {img_shape}")

    if len(img_shape) == 5:
        # (T, C, Z, Y, X) or similar
        t, c, z, y, x = img_shape
        if c != 1:
            raise ValueError("Image has multiple channels, please select TOM20 channel in advance")
    elif len(img_shape) == 4:
        # (C, Z, Y, X)
        c, z, y, x = img_shape
        if c != 1:
            raise ValueError("Image has multiple channels, please select TOM20 channel in advance")
        data = img[0]
    elif len(img_shape) == 3:
        data = img
    else:
        raise ValueError("Unsupported image dimension")

    # 保存为标准 Z-stack
    if output_path is None:
        filename = os.path.splitext(os.path.basename(input_path))[0]
        output_path = os.path.join(os.path.dirname(input_path), f"{filename}_zstack.tif")

    tiff.imwrite(output_path, data.astype(np.uint16))
    print(f"Saved cleaned Z-stack to: {output_path}")


# 示例用法（处理单个文件）
if __name__ == "__main__":
    input_file = '/Volumes/ExFAT/TFAM_fixed_JMB/TFAM_FFAA_slide8/tom20/TFAM_FFAA_slide8.tif'  # 修改为你的路径
    convert_hyperstack_to_zstack(input_file)
