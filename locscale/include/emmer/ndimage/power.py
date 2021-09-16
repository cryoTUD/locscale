import numpy as np
from scipy import fftpack

def psd2D_summed(images):
    psd2D_sum = None
    for image in images:
        if psd2D_sum is None:
            psd2D_sum = psd2D(image)
        else:
            psd2D_sum += psd2D(image)
    return np.array(psd2D_sum)

def psd2D(image):
    F1 = fftpack.fft2(image)
    F2 = fftpack.fftshift( F1 )
    psd2D = np.abs( F2 )**2
    return psd2D

def create_circular_mask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radius
    return mask.astype(float)

def create_circular_mask_3d(h, radius=None):
    center = int(h/2)
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center, h-center)
    Y, X, Z = np.ogrid[:h, :h, :h]
    dist_from_center = np.sqrt((X - center)**2 + (Y-center)**2 + (Z-center)**2)
    mask = dist_from_center <= radius
    return mask.astype(float)