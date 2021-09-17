import pandas
import gemmi
import numpy as np
import fcodes_fast
import emda.emda_methods as em
from emda.config import debug_mode


class Mtz:
    def __init__(self):
        self.filename = None
        self.h = None
        self.k = None
        self.l = None
        self.amplitudes = None
        self.phases = None
        self.f1d = None
        self.unit_cell = None
        self.f3d = None
        self.arr = None
        self.mapsize = None
        self.dataframe = None

    def _check_input(self):
        if self.filename is None:
            raise SystemExit("Filename is None")

    def read(self, mtzfilename=None):
        if self.filename is None:
            self.filename = mtzfilename
        self._check_input()
        try:
            mtz = gemmi.read_mtz_file(self.filename)
            self.unit_cell = np.zeros(6, dtype="float")
            self.unit_cell[0] = mtz.dataset(0).cell.a
            self.unit_cell[1] = mtz.dataset(0).cell.b
            self.unit_cell[2] = mtz.dataset(0).cell.c
            self.unit_cell[3:] = float(90)
            all_data = np.array(mtz, copy=False)
            dataframe = pandas.DataFrame(
                data=all_data, columns=mtz.column_labels())
            self.dataframe = dataframe
            for col in mtz.columns:
                if col.type == 'H' and col.label == 'H':
                    self.h = dataframe["H"].astype("int")
                if col.type == 'H' and col.label == 'K':
                    self.k = dataframe["K"].astype("int")
                if col.type == 'H' and col.label == 'L':
                    self.l = dataframe["L"].astype("int")
                if col.type == 'F':
                    self.amplitudes = dataframe[col.label]
                if col.type == 'P':
                    self.phase = dataframe[col.label]
        except FileNotFoundError as e:
            print(e)

    def to_f1d(self):
        if self.amplitudes is None:
            self.read()
        self.f1d = self.amplitudes * np.exp(np.pi * 1j * self.phase / 180.0)

    def to_f3d(self, mapsize=None):
        if self.f1d is None:
            self.to_f1d()
        if self.mapsize is None:
            self.mapsize = mapsize
        nx, ny, nz = self.mapsize
        self.f3d = fcodes_fast.mtz2_3d(self.h,
                                       self.k,
                                       self.l,
                                       self.f1d,
                                       nx, ny, nz,
                                       len(self.f1d))

    def to_map(self, mapsize=None):
        if self.mapsize is None:
            self.mapsize = mapsize
        if self.f3d is None:
            self.to_f3d()
        from numpy.fft import ifftshift, ifftn
        self.arr = np.real((ifftn(ifftshift(self.f3d))))


def mtz2map(mtzname, map_size):
    mtz = Mtz()
    mtz.read(mtzname)
    mtz.to_map(map_size)
    return mtz.arr, mtz.unit_cell


def mtz2f(mtzname):
    mtz = Mtz()
    mtz.read(mtzname)
    return mtz.to_f3d(), mtz.unit_cell


def _check_inputdata(mapdata):
    if isinstance(mapdata, np.ndarray):
        if mapdata.ndim == 3:
            if mapdata.dtype == np.complex128:
                return mapdata
            elif mapdata.dtype == np.float64:
                return np.fft.fftshift(np.fft.fftn(mapdata))
            else:
                print('data type: ', mapdata.dtype)
                raise SystemExit("wrong data type.")
    else:
        raise SystemExit("Not a numpy array")


def write_3d2mtz(unit_cell, mapdata, outfile="emda_map2mtz.mtz", resol=None):
    """ Writes 3D Numpy array into MTZ file.

    Arguments:
        Inputs:
            unit_cell: float, 1D array
                Unit cell params.
            mapdata: complex or float, 3D numpy array
                Map values to write.
            resol: float, optional
                Map will be output for this resolution.
                Default is up to Nyquist.

        Outputs: 
            outfile: string
            Output file name. Default is map2mtz.mtz.
    """
    from emda.core import restools

    mapdata = _check_inputdata(mapdata)
    nx, ny, nz = mapdata.shape
    assert nx == ny == nz
    nbin, res_arr, bin_idx = restools.get_resolution_array(unit_cell, mapdata)
    if resol is None:
        cbin = nx // 2
    else:
        dist = np.sqrt((res_arr - resol) ** 2)
        cbin = np.argmin(dist) + 1
        if cbin % 2 != 0:
            cbin += 1
        if nbin < cbin:
            cbin = nbin
    mtz = gemmi.Mtz()
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 1")
    mtz.cell.set(unit_cell[0], unit_cell[1], unit_cell[2], 90, 90, 90)
    mtz.add_dataset("HKL_base")
    for label in ["H", "K", "L"]:
        mtz.add_column(label, "H")
    mtz.add_column("Fout0", "F")
    mtz.add_column("Pout0", "P")
    # calling fortran function
    h, k, l, ampli, phase = fcodes_fast.prepare_hkl(
        mapdata, bin_idx, cbin, debug_mode, nx, ny, nz)
    print("hmin, hmax ", np.amin(h), np.amax(h))
    print("kmin, kmax ", np.amin(k), np.amax(k))
    print("lmin, lmax ", np.amin(l), np.amax(l))
    nrows = len(h)
    ncols = 5  # Change this according to what columns to write
    data = np.ndarray(shape=(nrows, ncols), dtype=object)
    data[:, 0] = -l.astype(int)
    data[:, 1] = -k.astype(int)
    data[:, 2] = -h.astype(int)
    data[:, 3] = ampli.astype(np.float32)
    data[:, 4] = phase.astype(np.float32)
    mtz.set_data(data)
    mtz.write_to_file(outfile)


if __name__ == "__main__":
    mtzname = '/Users/ranganaw/MRC/REFMAC/haemoglobin/EMD-3651/emda_test/map2mtz/map2mtz.mtz'
    map_size = [100, 100, 100]
    arr, uc = mtz2map(mtzname, map_size)
    em.write_mrc(arr, 'test2.mrc', uc)
