# Test code to get atomic correlation values by interpolation
# To use as a tool for local correlation example 2 in the tutorial
# 13.08.2021

import sys
import numpy as np
import gemmi
import emda.emda_methods as em


def interpolate_cc(data):
    from scipy.interpolate import RegularGridInterpolator

    linx = np.linspace(0, data.shape[0], data.shape[0])
    liny = np.linspace(0, data.shape[1], data.shape[1])
    linz = np.linspace(0, data.shape[2], data.shape[2])
    return RegularGridInterpolator(
        (linx, liny, linz), data, method="nearest")

    
def get_chains(modelname, cc_interpolator, pixsize):
    st = gemmi.read_structure(modelname)
    atomic_cc_list = []
    atom_name_list = []
    bval_list = []
    for chain in st[0]:
        for residue in chain:
            for atom in residue:
                if chain.name == 'C' and residue.name == 'EIC':
                    x = float(atom.pos.x) / pixsize
                    y = float(atom.pos.y) / pixsize
                    z = float(atom.pos.z) / pixsize
                    bval_list.append(float(atom.b_iso))
                    atomic_cc = cc_interpolator((z, y, x))
                    atomic_cc_list.append(atomic_cc)
                    atom_name_list.append(atom.name)
                    print(chain.name, residue.seqid.num, residue.name, atom.name, atomic_cc)
    return atomic_cc_list, atom_name_list, bval_list


def plot_cc(cclist, atom_name_list, biso_list):
    import numpy as np
    import matplotlib
    matplotlib.use(matplotlib.get_backend())
    from matplotlib import pyplot as plt
    from collections import deque
    from matplotlib.ticker import FormatStrFormatter

    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111)

    map_cc = np.roll(np.asarray(cclist[0], dtype='float'), 2)
    mapmodel_cc = np.roll(np.asarray(cclist[1], dtype='float'), 2)
    atomic_biso = np.roll(np.asarray(biso_list, dtype='float'), 2)
    atom_name_list_np = np.roll(np.asarray(atom_name_list), 2)
    items = deque(atom_name_list)
    cref = np.sqrt(map_cc)
    ax1.plot(mapmodel_cc, label='$CC_{map,model}$', marker='o', color='orange')
    ax1.plot(cref, label='$\sqrt{CC_{full}}$', marker='o', color='blue')

    ax2 = ax1.twinx()
    ax2.plot(atomic_biso, label="atomic Bs", marker='o', color='grey', alpha=0.5)
    
    x = np.arange(len(atom_name_list))
    y = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], dtype='float')
    yvals = np.array(["0.0", " ", "0.2", " ", "0.4", " ", "0.6", " ", "0.8", " ", "1.0"])
    print(items.rotate(2))
    plt.xticks(x, atom_name_list_np)
    ax1.set_xlabel("Atom label", fontweight='bold')
    ax1.set_ylabel("Correlation", fontweight='bold')
    ax1.set(yticks=y, yticklabels=yvals)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.legend(loc=3)
    ax1.grid(True, alpha=0.3, linestyle='--')
    # properties of ax2
    ax2.set_ylabel("Atomic B values", fontweight='bold')
    ax2.legend(loc=7)
    #
    plt.title("Atomic correlation values of linoleic acid", fontweight='bold')
    plt.savefig("atomic_cc.png", format="png", dpi=300)
    plt.show()


if __name__=="__main__":
    map_rcc, mapmodel_rcc, modelname = sys.argv[1:]

    rcc_map_list = [map_rcc, mapmodel_rcc]
    # read the correlation maps
    cc_list = []
    for i, imapcc in enumerate(rcc_map_list):
        print(imapcc)
        uc, rcc, orig = em.get_data(imapcc)
        pixsize = uc[0] / rcc.shape[0]
        cc_interpolator = interpolate_cc(data=rcc)
        if i == 0:
            cc, atom_name_list, biso_list = get_chains(modelname, cc_interpolator, pixsize)
            cc_list.append(cc)
        else:
            cc_list.append(get_chains(modelname, cc_interpolator, pixsize)[0])

    # plot various cc values
    plot_cc(cc_list, atom_name_list, biso_list)
