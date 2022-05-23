# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 10:06:12 2021
"""

import numpy as np
import pypdb
import gemmi
from locscale.include.emmer.pdb.pdb_tools import get_unit_cell_estimate, set_to_center_of_unit_cell
from locscale.include.emmer.pdb.pdb_utils import shift_coordinates
from locscale.include.emmer.ndimage.filter import low_pass_filter
from locscale.include.emmer.ndimage.map_utils import save_as_mrc

def detect_pdb_input(input_pdb):
    '''
    Function to detect the type of input the user has passed and return gemmi.Structure as output

    Parameters
    ----------
    input_pdb : str or gemmi.Structure
        Input can be either 
        a) string: either (i) pdb path, for ex: "path/to/pdb.pdb" or 
                          (ii)pdb_id "3j5p"
        b) gemmi.Structure()

    Returns
    -------
    pdb_structure : gemmi.Structure()
        Parsed input of type gemmi.Structure()
        Returns nothing if input cannot be parsed properly.

    '''
    from locscale.include.emmer.pdb.pdb_utils import get_gemmi_st_from_id
    
    if isinstance(input_pdb, str):
        if input_pdb.split(sep='.')[-1] in ['pdb','cif']:  # Then input is a file path
            pdb_structure = gemmi.read_structure(input_pdb)
            return pdb_structure
        elif ("/" not in input_pdb) or ("\\" not in input_pdb) or ("." not in input_pdb): # Then input is not file path but pdb_id
            pdb_structure = get_gemmi_st_from_id(input_pdb)
            return pdb_structure
        else:
            print("Input cannot be parsed. Please pass a pdb_id (as string) or pdb_path (as string)")
    elif isinstance(input_pdb, gemmi.Structure):
        pdb_structure = input_pdb.clone()
        return pdb_structure
    else:
        print("Unknown datatype for input. Please pass either a gemmi.Structure() or a\
              string (pointing to pdb_path or pdb_id")
    
              
def pdb2map(input_pdb=None, unitcell=None, size=None, apix=None, return_grid=False, verbose=False, 
            mdlidx=0,align_output=True,set_refmac_blur=True,set_unblur=True, blur=0):
    '''
    Cleaner function to convert a gemmi_structure to EM map. Make sure the input structure, or the pdb in the 
    path you input are correct. Common check include: 
        a) center of mass of the gemmi structure should be roughly in the center of unitcell (important)
        b) Remove waters 
        c) Check if atomic bfactors make sense
    
    Note: if you use a single value for apix then the program will assume same voxelsize in all dimensions

    Parameters
    ----------
    input_pdb : str or gemmi.Structure, required
        Input can be either 
        a) string: either (i) pdb path, for ex: "path/to/pdb.pdb" or 
                          (ii)pdb_id "3j5p"
        b) gemmi.Structure()

    unitcell : gemmi.UnitCell, 
        
    size : tuple, 
        Expected map shape
    apix : float, optional
        Expected voxelsize
    Either two of the three (unitcell/size/apix) allowed. In case only size is given as input, unitcell is taken 
    from gemmi structure
    return_grid : bool, optional
        
    verbose : TYPE, optional
        The default is False.
    mdlidx : int
        If gemmi.Structure has multiple models this index tells which model to select. Default is zero
    align_output : bool
        If selected, this transforms to output to align according to mrcfile convenctions 
        Transformation: flip axis: 2 then rotate in plane (2,0) by angle 90
        

    Returns
    -------
    if return_grid is set True: emmap (numpy.ndarray), grid (gemmi.FloatGrid)
    
    Else, only emmap is returned

    '''
    from scipy.ndimage import center_of_mass
    from locscale.include.emmer.ndimage.map_tools import sharpen_maps
    from locscale.include.emmer.ndimage.map_utils import convert_to_tuple
    pdb_structure = detect_pdb_input(input_pdb)
    
    if size is not None and apix is not None and unitcell is None:
        if verbose:
            print("Map size and pixelsize are given as inputs. Unitcell present in input structure will be ignored.")
        apix_gemmi = convert_to_tuple(apix, num_dims=3) # If apix is scalar, then the function returns a tuple of len 3
        unitcell_gemmi = gemmi.UnitCell(size[0]*apix_gemmi[0], size[1]*apix_gemmi[1], size[2]*apix_gemmi[2], 90, 90, 90)
        pdb_structure.cell = unitcell_gemmi
        size_gemmi = size
        
        
    
    if size is not None and unitcell is not None and apix is None:
        if verbose:
            print("Map size and unitcell are given as inputs. Unitcell present in input structure will be ignored.")
        unitcell_gemmi = unitcell
        pdb_structure.cell = unitcell_gemmi
        size_gemmi = size
        apix_gemmi = (unitcell_gemmi.a/size_gemmi[0],unitcell_gemmi.b/size_gemmi[1],unitcell_gemmi.c/size_gemmi[2])
    
    if unitcell is not None and apix is not None and size is None:
        if verbose:
            print("Pixel size and unitcell are given as inputs. Unitcell present in input structure will be ignored.")
        unitcell_gemmi = unitcell
        pdb_structure.cell = unitcell_gemmi
        apix_gemmi = convert_to_tuple(apix,num_dims=3)
        size_gemmi = [int(round(unitcell.a/apix_gemmi[0])), int(round(unitcell.b/apix_gemmi[1])), int(round(unitcell.c/apix_gemmi[2]))]
        
        
    if size is not None and apix is None and unitcell is None:  # in this case, program will use whatever unitcell is present in the gemmi structure
        size_gemmi = size
        apix_gemmi = (pdb_structure.cell.a/size_gemmi[0],pdb_structure.cell.b/size_gemmi[1],pdb_structure.cell.c/size_gemmi[2])

    reporter = {}
    reporter['pdb_struct_name'] =    pdb_structure.name
    reporter['unit_cell_exp'] =      pdb_structure.cell.parameters
    reporter['shape_exp'] = size_gemmi
    reporter['apix_exp'] = apix_gemmi
    reporter['com_pdb'] = np.array(pdb_structure[mdlidx].calculate_center_of_mass().tolist())
    if verbose:   # Make a check: 
        float_formatter = "{:.4f}".format
        np.set_printoptions(formatter={'float_kind':float_formatter})
        print("Simulating EM-map using: \t ", reporter['pdb_struct_name'], " \nwith the following properties: \n")
        print("Unit-cell (A*A*A): \t\t", reporter['unit_cell_exp'])
        print("Expected shape: \t\t", reporter['shape_exp'])
        print("Expected voxelsize (A*A*A): \t\t", reporter['apix_exp'])
        print("Center of mass: (A*A*A): \t\t",reporter['com_pdb'] )
        
        
    
    ## At this point, the values for size, pixelsize the unit cell are set 
    ## we should use the variable size_gemmi in density calculator

    if mdlidx > len(pdb_structure):  # Select model index from the pdb_structure
        print("selected model number not in pdb")
        return 0
    model = pdb_structure[mdlidx]
    dencalc = gemmi.DensityCalculatorE()
    dencalc.d_min = apix_gemmi[0]*2
    dencalc.rate = 1
    if blur > 0:
        dencalc.blur=blur
    if set_refmac_blur:
        dencalc.set_refmac_compatible_blur(model)
        inv_d2 = dencalc.blur
        if verbose:
            print("Setting a blur of {:.2f}".format(dencalc.blur))
            
    dencalc.set_grid_cell_and_spacegroup(pdb_structure)
    try:
        dencalc.grid.set_size(*size_gemmi)
        dencalc.add_model_density_to_grid(model)
    except:
        dencalc.put_model_density_on_grid(model)

    emmap = np.array(dencalc.grid,copy=False)
    
    if set_refmac_blur:
        if verbose:
            print("Applying a unblur for the sampled density equal to: {:.2f}".format(-inv_d2))
        emmap = sharpen_maps(emmap, apix=apix_gemmi[0], global_bfactor=-inv_d2)
        
    
    if align_output:
        from scipy.ndimage import rotate
        emmap_flipped = np.flip(emmap,axis=2)
        emmap_rotated = rotate(emmap_flipped, angle=90, axes=(2,0))
        emmap = emmap_rotated
    
    reporter['shape_final'] = emmap.shape
    reporter['apix_final'] = np.array(dencalc.grid.spacing)
    reporter['com_map'] = np.array(center_of_mass(abs(emmap)))*reporter['apix_final']
    if verbose: ## Check output if it matches expectation
        print("\nMap simulated! Final parameters:")
        print("Emmap shape \t\t",reporter['shape_final'])
        print("Grid spacing (A*A*A): \t\t", reporter['apix_final'])
        print("Center of mass of Emmap: (A*A*A)\t", reporter['com_map'])
    
    
    if return_grid:
        return emmap, dencalc.grid
    
    else:
        return emmap
    
