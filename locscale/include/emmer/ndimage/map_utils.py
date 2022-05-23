import numpy as np

def load_map(map_path, return_apix = True):
    import mrcfile
    from locscale.include.emmer.ndimage.map_utils import average_voxel_size
    emmap = mrcfile.open(map_path).data
    apix = average_voxel_size(mrcfile.open(map_path).voxel_size)
    
    if return_apix:
        return emmap, apix
    else:
        return emmap
def parse_input(input_map, allow_any_dims=True):
    '''
    Function to detect type of input and return a emmap numpy array

    Parameters
    ----------
    input_map : str or numpy array or Mrc object
    string type input should be path/to/emmap.mrc    

    Returns
    -------
    emmap : numpy.ndarray

    '''
    import os
    import mrcfile
    if isinstance(input_map, np.ndarray):
        if not allow_any_dims:
            if len(input_map.shape) == 3:
                return input_map
            else:
                print("You have not input a 3-D numpy array, which cannot be a EM-map")
                return None
        else:
            return input_map
            
    elif isinstance(input_map, str):
        if os.path.exists(input_map):
            emmap = mrcfile.open(input_map).data
            return emmap
        else:
            print("You have not entered a proper path, or the requested file does not exist!")
            return None
    
    
            
def get_all_voxels_inside_mask(mask_input, mask_threshold=1):
    mask = parse_input(mask_input)
    all_inside_mask = np.asarray(np.where(mask>=mask_threshold)).T.tolist()
    return all_inside_mask
    
    
def read_gemmi_map(map_path, return_grid=False):
    '''
    Function to read a map file and return a numpy.ndarray using gemmi.read_ccp4() map

    Parameters
    ----------
    map_path : str
        path/to/emmap.mrc

    Returns
    -------
    emmap : numpy.ndarray
        
    '''
    import gemmi
    
    gemmi_ccp4Map = gemmi.read_ccp4_map(map_path)
    emmap = np.array(gemmi_ccp4Map.grid, copy=False)
    
    if return_grid:
        return emmap, gemmi_ccp4Map.grid
    else:
        return emmap
    
    
def save_as_mrc(map_data,output_filename, apix=None,origin=None,verbose=False, header=None):
    '''
    Function to save a numpy array containing volume, as a MRC file with proper header

    Parameters
    ----------
    map_data : numpy.ndarray
        Volume data showing the intensities of the EM Map at different points

    apix : float or any iterable
        In case voxelsize in x,y,z are all equal you can also just pass one parameter. 
    output_filename : str
        Path to save the MRC file. Example: 'path/to/map.mrc'
    origin: float or any iterable, optional
        In case origin index in x,y,z are all equal you can also just pass one parameter. 

    Returns
    -------
    Saves MRC .

    '''
    import mrcfile

    with mrcfile.new(output_filename,overwrite=True) as mrc:
        mrc.set_data(np.float32(map_data))
        
        if header is not None:
            mrc.set_extended_header(header)
        
        else:
            if apix is not None:
                #apix_list = [apix['x'], apix['y'], apix['z']]
                ## apix can be either a float or a list. If it's a single number, then the function convert_to_tuple will use it three times
                apix_tuple = convert_to_tuple(apix, num_dims=3)
                rec_array_apix = np.rec.array(apix_tuple, dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
                mrc.voxel_size = rec_array_apix
            else:
                print("Please pass a voxelsize value either as a float or an iterable")
                return 0
            
            if origin is not None:    
                origin_tuple = convert_to_tuple(origin,num_dims=3)
            else:
                origin_tuple = convert_to_tuple(input_variable=0,num_dims=3)
            rec_array_origin = np.rec.array(origin_tuple, dtype=[('x','<f4'),('y','<f4'),('z','<f4')])
            mrc.header.origin = origin_tuple
            
        if verbose:
            print("Saving as MRC file format with following properties: ")
            print("Voxel size", mrc.voxel_size)
            print("Origin", mrc.header.origin)
        
    mrc.close()


def compare_gemmi_grids(grid_1, grid_2):
    '''
    Function to test similarity of two gemmi grids. Test include: 
        (a) Axis order
        (b) Voxelsize
        (c) UnitCell
        (d) Shape
        

    Parameters
    ----------
    grid_1 : gemmi.FloatGrid
        
    grid_2 : gemmi.FloatGrid

    Returns
    -------
    report : pandas.DataFrame
    

    '''
    import pandas as pd
    report = pd.DataFrame()
    report['axis_order'] = [grid_1.axis_order.name,grid_2.axis_order.name]
    report['spacing'] = [grid_1.spacing,grid_2.spacing]
    report['unitcell'] = [grid_1.unit_cell,grid_2.unit_cell]
    report['shape'] = [grid_1.shape, grid_2.shape]
    
    report = report.T
    report['final'] = report[0] == report[1] 
    # if report['final'].all():
    #     print("The two input grids are same")
    # else:
    #     print("Two input grids are not the same")
    #     print(report['final'])
    return report

def ZYX_to_XYZ(emmap):
    '''
    Function to convert a ZYX numpy array to XYZ numpy array

    Parameters
    ----------
    emmap : numpy.ndarray
        Volume data showing the intensities of the EM Map at different points

    Returns
    -------
    emmap : numpy.ndarray
        
    '''
    from scipy.ndimage import rotate
    emmap_flipped = np.flip(emmap,axis=2)
    emmap_rotated = rotate(emmap_flipped, angle=90, axes=(2,0))

    return emmap_rotated    
def extract_window(im, center, size):
    '''
    Extract a square window at a given location. 
    The center position of the window should be provided.

    Parameters
    ----------
    im : numpy.ndarray
        3D numpy array
    center : tuple, or list, or numpy.array (size=3)
        Position of the center of the window
    size : int, even
        Total window size (edge to edge) as an even number
        (In future could be modified to include different sized window 
        in different directions)
        

    Returns
    -------
    window : numpy.ndarray
        3D numpy array of shape (size x size x size)

    '''
    z,y,x = center
    window = im[z-size//2:z+size//2, y-size//2:y+size//2, x-size//2:x+size//2]
    return window

def binarizeMap(emmap, threshold):
    binary_map = (emmap>=threshold).astype(np.int_)
    return binary_map

def average_voxel_size(voxel_size_record):
    apix_x = voxel_size_record.x
    apix_y = voxel_size_record.y
    apix_z = voxel_size_record.z
    
    average_apix = (apix_x+apix_y+apix_z)/3
    
    return average_apix

def compute_FDR_confidenceMap_easy(em_map, apix, window_size, fdr=1, lowPassFilter_resolution=None,remove_temp_files=True, folder = None):
    from locscale.include.confidenceMapUtil.confidenceMapMain import calculateConfidenceMap
    import os, shutil, time
    
    if folder is None:
        current_cwd = os.getcwd()
    else:
        current_cwd = folder
    timestamp =  str(time.time())
    temp_dir = current_cwd + '/fdr_output_temp_'+timestamp
    os.mkdir(temp_dir)
    os.chdir(temp_dir)
    confidenceMap,locFiltMap,locScaleMap,binMap,maskedMap = calculateConfidenceMap(
        em_map=em_map,apix=apix,noiseBox=None,testProc=None,ecdf=None,
        lowPassFilter_resolution=lowPassFilter_resolution,method=None, 
        window_size=window_size,windowSizeLocScale=None, locResMap=None,
        meanMap=None,varMap=None,fdr=fdr,modelMap=None,stepSize=None,mpi=None)
    
    fdr_threshold = np.min(maskedMap[np.nonzero(maskedMap)])
    
    os.chdir(current_cwd)
    if remove_temp_files:
        print("Clearing temporary files")
        shutil.rmtree(temp_dir)
    return confidenceMap, fdr_threshold
    
def get_sphere(radius):
    '''
    Function to return a window, with a spherical mask. Size of the window defined by the radius.
    Parameters
    ----------
    radius : int
        Radius of sphere, in pixels

    Returns
    -------
    sphere : numpy.ndarray
        Shape: (rad+1) x (rad+1)

    '''
    z,y,x = np.ogrid[-radius: radius+1, -radius: radius+1, -radius: radius+1]
    sphere = (x**2+y**2+z**2 <= radius**2).astype(int)
    return sphere


def pad_or_crop_image(im, pad_factor=None, pad_value = None, crop_image=False):
    """Returns the original image cropped or padded by pad_factor and pad_value; pad_factor being a fraction/multiple of original image size.
       Default behaviour is zero padding.
    """
    if np.any(pad_factor == None):
        return im
    else:
        pad_factor = np.round(np.multiply(pad_factor,np.array(im.shape))).astype('int')

        if pad_value == None:
            pad_value = 0

        if len(im.shape) == 2:       
            if (pad_factor[0] <= im.shape[0] or pad_factor[1] <= im.shape[1]):
                crop_image = True    
            
            if crop_image:
                crop_im = im[im.shape[0]//2-pad_factor[0]//2:im.shape[0]//2+pad_factor[0]//2+pad_factor[0]%2, :]
                crop_im = crop_im[:, im.shape[1]//2-pad_factor[1]//2:im.shape[1]//2+pad_factor[1]//2+pad_factor[1]%2]
                return crop_im
            else:
                pad_im = np.pad(im, ((pad_factor[0]//2-im.shape[0]//2, pad_factor[0]//2-im.shape[0]//2+pad_factor[0]%2), (0,0)), 'constant', constant_values=(pad_value,))
                pad_im = np.pad(pad_im, ((0,0),(pad_factor[1]//2-im.shape[1]//2, pad_factor[1]//2-im.shape[1]//2+pad_factor[1]%2 )), 'constant', constant_values=(pad_value,))
                return pad_im         
            
        elif len(im.shape) == 3:
            if (pad_factor[0] <= im.shape[0] or pad_factor[1] <= im.shape[1] or pad_factor[2] <= im.shape[2]):
                crop_image = True

            if crop_image:
                crop_im = im[im.shape[0]//2-pad_factor[0]//2:im.shape[0]//2+pad_factor[0]//2+pad_factor[0]%2, :, :]
                crop_im = crop_im[:, im.shape[1]//2-pad_factor[1]//2:im.shape[1]//2+pad_factor[1]//2+pad_factor[1]%2, :]
                crop_im = crop_im[:, :, im.shape[2]//2-pad_factor[2]//2:im.shape[2]//2+pad_factor[2]//2+pad_factor[2]%2]
                return crop_im

            else:
                pad_im = np.pad(im, ((pad_factor[0]//2-im.shape[0]//2, pad_factor[0]//2-im.shape[0]//2+pad_factor[0]%2), (0,0), (0,0) ), 'constant', constant_values=(pad_value,))
                pad_im = np.pad(pad_im, ((0,0), (pad_factor[1]//2-im.shape[1]//2, pad_factor[1]//2-im.shape[1]//2+pad_factor[1]%2 ), (0,0)), 'constant', constant_values=(pad_value,))
                pad_im = np.pad(pad_im, ((0,0), (0,0), (pad_factor[2]//2-im.shape[2]//2, pad_factor[2]//2-im.shape[2]//2+pad_factor[2]%2)), 'constant', constant_values=(pad_value,))
                return pad_im

def resample_image(im, imsize_new=None, apix=1.0, apix_new=None):
    """Returns a real image or volume resampled by cropping/padding its Fourier Transform
    """
    import numpy as np
    
    imsize = np.array(im.shape)
    if np.any(imsize_new == None) and apix_new == None:
        imsize_new = im.shape
        apix_new = apix
        pad_factor = imsize_new[0]/imsize[0]
        pad_factor = np.round(np.array(tuple([pad_factor*i for i in im.shape]))).astype('int')
    elif apix_new != None:
        imsize_new = np.round(imsize * apix / apix_new).astype('int')
        pad_factor = imsize_new/imsize
    elif imsize_new != None:
        imsize_new = np.array(imsize_new)
        pad_factor = imsize_new/imsize
    
    
    ft = np.fft.fftn(im)
    ft = np.fft.fftshift(ft)
    
    ft = pad_or_crop_image(ft, pad_factor)
    
    ft = np.fft.ifftshift(ft)
    
    real_image =np.fft.ifftn(ft).real
    return real_image

def resample_map(emmap, emmap_size_new=None, apix=None, apix_new=None, order=1):
    '''
    Function to resample an emmap in real space using linear interpolation 

    Parameters
    ----------
    emmap : numpy.ndimage
        
    emmap_size_new : tuple 
        
    apix : float
        
    apix_new : float
        

    Returns
    -------
    resampled_emmap

    '''
    from scipy.ndimage import zoom
    if emmap_size_new is None:
        if apix is not None and apix_new is not None:
            resample_factor = apix/apix_new
        else:
            raise UserWarning("Provide either (1) current pixel size and new pixel size or (2) new emmap size")
    
    else:
        try:
            resample_factor = emmap_size_new[0] / emmap.shape[0]
        except:
            raise UserWarning("Please provide proper input: emmap_size_new must be a tuple")
    
    resampled_image = zoom(emmap, resample_factor, order=order)
    
    return resampled_image


def measure_mask_parameters(mask_path=None, mask=None,apix=None,edge_threshold=0.99,protein_density=1.35,average_atomic_weight=13.14,verbose=True,detailed_report=False):
    import mrcfile
    import numpy as np
    from scipy.constants import Avogadro
    '''
    Function to calculated parameters of a EM Mask map

    Parameters
    ----------
    mask_path : string 
        Path to mask file
    edge_threshold : float 
        The threshold to strictly binarize the FDR map at the edges
    protein_density : float, optional
        Average protein density to calculate number of atoms. The default is 1.35.
    average_atomic_weight : float, optional
        Atomic weight of an "average atom present in protein". 
        Found using 54% carbon, 20% oxygen and 16% nitrogen. The default is 12.066.
    verbose : bool, optional
        Print statistics if True. The default is True.

    Returns
    -------
    num_atoms : int
        Estimated number of atoms based on mask volume, protein density and average atomic weight
    

    '''
    if mask_path is not None:
        mask_mrc = mrcfile.open(mask_path)
        mask = mask_mrc.data
        voxelsize = mask_mrc.voxel_size.x
    elif mask is not None and apix is not None:
        mask = mask
        voxelsize = apix
    elif mask_path is None or (mask is None and apix is None):
        print("Input error: Provide atleast mask path, or (mask and apix)")
        return None
    
    ang_to_cm = 1e-8
    
    mask = (mask>=edge_threshold).astype(np.int_)
    
    mask_vol = mask.sum()*(voxelsize*ang_to_cm)**3
    mask_vol_A3 = mask.sum()*voxelsize**3
    #print("\n Volume of the mask generated is: "+str(mask.sum())+" A$^3$ \n")
    # Calculate number of atoms
    protein_mass = protein_density * mask_vol
    num_moles = protein_mass / average_atomic_weight
    num_atoms = int((num_moles * Avogadro).round())
    maskshape = mask.shape
    mask_dims = [maskshape[0]*voxelsize,maskshape[1]*voxelsize,maskshape[2]*voxelsize]
    
    
    if verbose:
        print("Mask parameters calculated are: \n"+
              "Mask sum voxels: "+str(round(mask.sum(),3))+"\n"+
              "Mask volume: "+str(round(mask_vol_A3,3))+" A^3 \n"+
              "Protein mass: "+str(round(1e21*protein_mass))+" zg\n"+
              "Num atoms: "+str(num_atoms)+"\n") 
        
    if not detailed_report:
        return num_atoms,mask_dims
    else:
        return mask_vol_A3, protein_mass, num_atoms, mask_dims,maskshape
    

def get_nyquist_limit(xdata):
    dx = xdata[1]-xdata[0]
    return 1/(dx*2)


def fit_series(series,xmin,xmax,num):
    from scipy.interpolate import interp1d
        
    xdata = series[0]
    ydata = series[1]
    
 #   print("Fitting series. Initial frequency range:  "+str(tuple([xdata[0],xdata[-1]]))+"\t Final frequency range: "+str(tuple([xmin,xmax])))
    f = interp1d(xdata,ydata,fill_value='extrapolate')
    new_xdata = np.linspace(xmin,xmax,num)
    new_ydata = f(new_xdata)
    
    '''
    
    if xdata[0] > xmin:
        print("fitting left half")
        f_left = interp1d(xdata[:3],ydata[:3],fill_value='extrapolate')
        dx_left = (xdata[2] - xdata[0])/2
        left_xnew = np.arange(xmin,xdata[0],dx_left)
        left_ynew = f(left_xnew)
        xdata = np.concatenate((left_xnew,xdata))
        ydata = np.concatenate((left_ynew,ydata))
    if xdata[-1] < xmax:
        print("fitting right half")
        f_right = interp1d(xdata[-3:],ydata[-3:],fill_value='extrapolate')
        dx_right = (xdata[-3] - xdata[-1])/2
        right_xnew = np.arange(xdata[-1],xmax,dx_right)
        right_ynew = f(right_xnew)
        xdata = np.concatenate((xdata,right_xnew))
        ydata = np.concatenate((ydata,right_ynew))
    
    '''
    return [new_xdata,new_ydata]

def moving_average(array,window=5):
    return np.convolve(array,np.ones(window), 'same')/window

def normalise(x):
    if isinstance(x,np.ndarray):
        return x/x.max()
    else:
        x = np.array(x)
        return x/x.max()
    
def convert_pdb_to_mrc_position(pdb_position, apix):
    '''
    Convert the real units of positions into indices for the emmap. 
    Note: returns in (Z,Y,X) format
    

    Parameters
    ----------
    pdb_position : list
        list of xyz positions (Angstorm)
    apix : float
        Pixel size 

    Returns
    -------
    mrc_position : list
        List of ZYX positions (index positions)

    '''
    mrc_position = []
    
    for pos in pdb_position:
        [x,y,z] = pos
        int_x, int_y, int_z = int(round(x/apix)), int(round(y/apix)), int(round(z/apix))
        mrc_position.append([int_z, int_y, int_x])
        
    return mrc_position

def convert_mrc_to_pdb_position(mrc_position_list, apix):
    '''
    Convert the real units of positions into indices for the emmap. 
    Note: returns in (Z,Y,X) format
    
    Parameters
    ----------
    mrc_position_list : list
        list of xyz positions (Angstorm)
    apix : float
        Pixel size 

    Returns
    -------
    pdb_position_list : list
        List of XYZ positions (index positions)

    '''
    pdb_position_list = []
    
    for pos in mrc_position_list:
        [nz,ny,nx] = pos
        z, y, x  = nz*apix, ny*apix, nx*apix
        pdb_position_list.append([x, y, z])
        
    return pdb_position_list

def dilate_mask(mask, radius, iterations=1):
    '''
    Dilate mask with spherical structures

    Parameters
    ----------
    mask : numpy.ndarray
        Skeleton structure of a set of atoms
    radius : int
        Cutoff radius used for binary dilation
    iterations : int, optional
        Number of iterations for binary dilation The default is 1.

    Returns
    -------
    dilated : numpy.ndarray
        Dilated mask

    '''
    from scipy.ndimage import binary_dilation
    from locscale.include.emmer.ndimage.map_utils import get_sphere
    
    dilated = binary_dilation(mask, structure=get_sphere(radius), iterations=iterations).astype(int)
        
    return dilated

def convert_to_tuple(input_variable, num_dims=3):
    '''
    Convert any variable, or iterable into a tuple. If a scalar is input then a tuple is generated with same variable
    based on number of dimensions mentioned in num_dims

    Parameters
    ----------
    input_variable : any
        scalar, or any iterable
    num_dims : int, optional
        Length of tuple. The default is 3.
        
    Returns
    -------
    output_tuple : tuple

    '''
    
    if hasattr(input_variable, '__iter__'):
        if len(input_variable) == num_dims:
            output_tuple = tuple(input_variable)
            return output_tuple
        else:
            print("Input variable dimension {} doesn't match expected output dimension {}".format(len(input_variable), num_dims))
    else:
        output_list = [input_variable for temporary_index in range(num_dims)]
        output_tuple = tuple(output_list)
        return output_tuple
    
def convert_to_tuple_2(input_variable, num_dims=3):
    '''
    Convert any variable, or iterable into a tuple. If a scalar is input then a tuple is generated with same variable
    based on number of dimensions mentioned in num_dims

    Parameters
    ----------
    input_variable : any
        scalar, or any iterable
    num_dims : int, optional
        Length of tuple. The default is 3.
        
    Returns
    -------
    output_tuple : tuple

    '''
    
    if hasattr(input_variable, '__iter__'):
        if isinstance(input_variable, np.recarray):
            input_variable = [input_variable['x'], input_variable['y'], input_variable['z']]
        
        if len(input_variable) == num_dims:
            output_tuple = tuple(input_variable)
            return output_tuple
        else:
            print("Input variable dimension {} doesn't match expected output dimension {}".format(len(input_variable), num_dims))
    else:
        output_list = [input_variable for temporary_index in range(num_dims)]
        output_tuple = tuple(output_list)
        return output_tuple  

def get_model_mask(input_pdb, mask_shape, apix, threshold_factor=5, bfactor=100, smoothen=5, tight=0.9):
    '''
    Simulates a map from an atomic model at a high bfactor and binarises it to result in a mask

    Parameters
    ----------
    input_pdb : TYPE
        DESCRIPTION.
    mask_threshold : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    mask : numpy.ndarray

    '''
    
    import gemmi
    from locscale.include.emmer.pdb.pdb_to_map import detect_pdb_input, pdb2map
    from locscale.include.emmer.pdb.pdb_utils import set_atomic_bfactors
    from scipy.signal import fftconvolve
    
    gemmi_st = detect_pdb_input(input_pdb)
    
    gemmi_st_bfactor = set_atomic_bfactors(input_gemmi_st=gemmi_st, b_iso=bfactor)
    
    simmap = pdb2map(input_pdb=gemmi_st_bfactor, apix=apix, size=mask_shape, set_refmac_blur=True, verbose=True)
    
    simmap_array = simmap.flatten()
    nonzero_array = simmap_array[simmap_array>0]
    threshold = np.percentile(nonzero_array, threshold_factor)
    binarised_map = (simmap>=threshold).astype(np.int_)
    
    kernel = np.ones((smoothen,smoothen, smoothen))
    smoothened_map = fftconvolve(in1=binarised_map, in2=kernel, mode="same")
    smoothened_map = smoothened_map/smoothened_map.max()
    
    binarised_map = (smoothened_map>=tight).astype(np.int_)
    
    return binarised_map
    
    
