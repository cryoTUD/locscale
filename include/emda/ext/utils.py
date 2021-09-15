# various functions
import numpy as np
import fcodes_fast

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def regularize_fsc_weights(fsc):
    for i, elm in enumerate(fsc):
        if elm <= 1e-5:
            elm = elm + 1e-4 
            fsc[i] = elm
    return fsc

def get_avg_fsc(binfsc, bincounts):
    fsc_filtered = filter_fsc(bin_fsc=binfsc, thresh=0.)
    fsc_avg = np.average(a=fsc_filtered[np.nonzero(fsc_filtered)], 
                         weights=bincounts[np.nonzero(fsc_filtered)])
    return fsc_avg

def filter_fsc(bin_fsc, thresh=0.1):
    bin_fsc_new = np.zeros(bin_fsc.shape, 'float')
    for i, ifsc in enumerate(bin_fsc):
        if ifsc >= thresh:
            bin_fsc_new[i] = ifsc
        else:
            if i > 1:
                break
    return bin_fsc_new

def halve_the_dim(arr1):
    nx, ny, nz = arr1.shape
    dx = int(nx / 4)
    dy = int(ny / 4)
    dz = int(nz / 4)
    return arr1[dx:dx+nx//2, dy:dy+ny//2, dz:dz+nz//2]

def get_ibin(bin_fsc, cutoff):
    # search from rear end
    ibin = 0
    for i, ifsc in reversed(list(enumerate(bin_fsc))):
        if ifsc > cutoff:
            ibin = i
            if ibin % 2 != 0:
                ibin = ibin - 1
            break
    return ibin


def determine_ibin(bin_fsc, cutoff=0.20):
    bin_fsc = filter_fsc(bin_fsc)
    bin_fsc = bin_fsc[np.nonzero(bin_fsc)]
    cutoff = np.average(bin_fsc) * 0.1
    ibin = get_ibin(bin_fsc, cutoff)        
    """ i = 0
    while ibin < 5:
        cutoff -= 0.01
        ibin = get_ibin(bin_fsc, max([cutoff, 0.1]))
        i += 1
        if i > 200:
            print("Fit starting configurations are too far.")
            raise SystemExit() """
    if ibin == 0:
        print("Fit starting configurations are too far.")
        raise SystemExit()
    return ibin


def fsc_between_static_and_transfomed_map(maps, bin_idx, rm, t, nbin):
    import emda.core as core

    nx, ny, nz = maps[0].shape
    maps2send = np.stack((maps[1], maps[2]), axis = -1)
    #frs = fcodes_fast.trilinear2(maps2send,bin_idx,rm,nbin,0,2,nx, ny, nz)
    frs = fcodes_fast.tricubic(rm=rm,
                                f=maps2send,
                                mode=0,
                                nc=2,
                                nx=nx, ny=ny, nz=nz) 
    st, _, _, _ = fcodes_fast.get_st(nx, ny, nz, t)
    f1f2_fsc, _, bin_count = core.fsc.anytwomaps_fsc_covariance(
        maps[0], frs[:, :, :, 0] * st, bin_idx, nbin)
    fsc_avg = get_avg_fsc(binfsc=f1f2_fsc, bincounts=bin_count)
    return [f1f2_fsc, frs[:, :, :, 0] * st, frs[:, :, :, 1] * st, fsc_avg]