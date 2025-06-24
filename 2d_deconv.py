import numpy as np
from skimage.io import imread, imsave
import RedLionfishDeconv as rl
from tqdm import tqdm

img_path = r'F:\atp\deconvolution\data\1.TIF'
psf_path = r'F:\atp\PSF BW.tif'
output_path = r'F:\atp\deconvolution\deconv_1_slice18_iters.tif'

image = imread(img_path)
psf = imread(psf_path)

# 取第18和19个slice（索引17和18）
img_3d = image[17:19]
psf_3d = psf[17:19]

results = []

for n_iter in tqdm(range(1, 101), desc='Deconvolving'):
    deconv = rl.doRLDeconvolutionFromNpArrays(
        img_3d, psf_3d, niter=n_iter, method='gpu', resAsUint8=False
    )
    if deconv is not None:
        # 只取第0层（即原本的第18层）
        results.append(deconv[0].astype(np.float32))
    else:
        print(f'迭代{n_iter}时反卷积返回None，跳过')
        results.append(np.zeros_like(img_3d[0], dtype=np.float32))

results_stack = np.stack(results, axis=0)
imsave(output_path, results_stack)
print(f'保存完成: {output_path}')