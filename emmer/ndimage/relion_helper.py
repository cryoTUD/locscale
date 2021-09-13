import mrcfile
from scipy.ndimage import rotate
import numpy as np

def load_particle_images_from_star(names, angles=None, prefix=""):
    images = []
    for i,name in enumerate(names):
        number,path = name.split('@')
        number = int(number)
        with mrcfile.open(prefix+path) as mrc:
            image = mrc.data[number-1]
            image = image - image.mean()
            image /= image.std()
            #image -= gaussian_filter(image, 30)
            if angles is not None:
                angle = angles[i]
                image = rotate(image, -angle+90, reshape=False)
            images.append(image)
    return np.array(images)