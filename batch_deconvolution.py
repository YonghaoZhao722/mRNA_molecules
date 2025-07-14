import os
from skimage.io import imread, imsave
import numpy as np
import RedLionfishDeconv as rl
import tqdm

data_dir = r'../mitochondria_FITC/Y333 ATP6 TIM50/FITC'
psf_path = r'../mitochondria_FITC/PSF BW.tif'
output_dir = r'Y333 ATP6 TIM50/deconvolved_30'
os.makedirs(output_dir, exist_ok=True)

psf = imread(psf_path)

iterations = 30

for fname in tqdm.tqdm(os.listdir(data_dir), desc='Processing', total=len(os.listdir(data_dir))):
    if fname.lower().endswith('.tif') and not fname.startswith('deconv_'):
        img_path = os.path.join(data_dir, fname)
        try:
            image = imread(img_path)
            deconvolved = rl.doRLDeconvolutionFromNpArrays(
                image, psf, niter=iterations, method='gpu', resAsUint8=False,
            )
            out_path = os.path.join(output_dir, f'deconv_{fname}')
            imsave(out_path, deconvolved.astype(np.float32))
            # print(f'Saved: {out_path}')
        except Exception as e:
            print(f'Error processing {img_path}: {e}')
print('All done!')