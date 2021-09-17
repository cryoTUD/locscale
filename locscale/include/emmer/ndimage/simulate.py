import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage


def make_2d_letter(text, image_size=144):
    fig,ax = plt.subplots(figsize=(2,2))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis('off')

    fig.text(0.23,0.26,text, fontsize=100)
    fig.canvas.draw()

    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    plt.close()
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,)).mean(2)
    data -= data.min(); data /= data.max()
    data = data.astype(np.bool)[::-1]
    data_inv = (~data).astype(float)
    data_zoom = scipy.ndimage.zoom(data_inv, image_size/data_inv.shape[0])
    return data_zoom 
