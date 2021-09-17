def fakehalf(f_map, knl=3):
    """Average neighbor Fourier Coefficients to reconstruct other half."""
    import numpy as np
    from emda import core
    import scipy.signal

    kernel = core.restools.create_soft_edged_kernel_pxl(knl)
    kernel[kernel.shape[0] // 2, kernel.shape[1] // 2, kernel.shape[2] // 2] = 0
    kernel = kernel / np.sum(kernel)
    fakehalf = scipy.signal.fftconvolve(f_map, kernel, "same")
    return fakehalf


def get_index(arr):
    """Finds the array index which has the lowest value.

    Search begins from rear end and continue forward.

    """
    f1 = False
    for i, _ in enumerate(arr):
        if i == 0:
            val_cur = arr[-1]
            val_pre = arr[-1]
        else:
            j = -i - 1
            val_cur = arr[j]
            if val_cur >= val_pre:
                if not f1:
                    val_pre = val_cur
                if f1:
                    index = len(arr) + j + 1
                    break
            elif val_cur < val_pre:
                f1 = True
                val_pre = val_cur
    return index