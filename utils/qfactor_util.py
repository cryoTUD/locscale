#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 13:17:09 2021

@author: alok
"""
import numpy as np

class EMDATA():
    def __init__(self,map_data):
        self.data = map_data
    
    def extract_window_from_index(self,center,size):
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
        window = self.data[z-size//2:z+size//2, y-size//2:y+size//2, x-size//2:x+size//2]
        return window
    
    def extract_window_from_position(self,position, distance, apix, round_up='odd'):
        '''
        Extract window from a position in Angstorm

        Parameters
        ----------
        position : TYPE
            DESCRIPTION.
        distance : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        from emmer.ndimage.map_utils import convert_pdb_to_mrc_position
        from utils.general import round_up_to_odd, round_up_to_even
        
        
        mrc_position = convert_pdb_to_mrc_position([position], apix)[0]
        
        if round_up == 'odd':
            size = round_up_to_odd(distance / apix)
        elif round_up == 'even':
            size = round_up_to_even(distance / apix)
        
        window = self.extract_window_from_index(mrc_position, size)
        
        return window

    
def compute_radial_intensity(vol, center=None, return_indices=False):
    '''
    Computes the radial profile of a given volume

    Parameters
    ----------
    vol : numpy.ndarray
        Input array
    center : list, optional
        DESCRIPTION. The default is [0,0,0].
    return_indices : bool, optional
        

    Returns
    -------
    radial_profile : numpy.ndarray (1D)
        Radial profile
        

    '''
    dim = vol.shape
    m = np.mod(vol.shape,2)

    if center is None:

        z, y, x = np.indices(vol.shape)
        center = tuple((a - 1) / 2.0 for a in vol.shape[::-1])
        radii = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
        radii = radii.astype(np.int)
    
    radial_profile = np.bincount(radii.ravel(), vol.ravel()) / np.bincount(radii.ravel())
    radial_profile[0:int(vol.shape[0]/2)]
    if not return_indices:
        return radial_profile
    else:
        return radial_profile, radii
