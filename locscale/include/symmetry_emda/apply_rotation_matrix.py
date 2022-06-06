## Script to perform trilinear interpolation operation on a 3D FT image
import numpy as np
def get_spherical_mask(shape, radius_pixels):
    if isinstance(shape, int):
        shape = (shape, shape, shape)

    n = shape[0]    
    z,y,x = np.ogrid[-n//2:n//2,-n//2:n//2,-n//2:n//2]
    mask = (x**2+y**2+z**2 <= radius_pixels**2).astype(np.int)
    return mask


def transform_origin_to_center_of_image(coordinates, center, shape):
    """
    Transform a set of 3D numpy indices to XYZ coordinates, with the origin at the center of the image
    """
    import numpy as np
    coordinates = np.array(coordinates)
    center = np.array(center)
    transformed_coordinates = coordinates - center  ## still in ZYX
    # convert to XYZ
    xyz_transformed_coordinates = np.flip(transformed_coordinates,axis=1)
    # Invert y axis
    invert_matrix = np.ones(transformed_coordinates.shape)
    invert_matrix[:,1] = -1
    invert_y_axis = xyz_transformed_coordinates * invert_matrix
    
    return np.clip(invert_y_axis, 0, shape[0]//2-1)

def transform_coordinates_to_numpy_coordinates(coordinates, center, shape):
    """
    Transform a set of 3D numpy indices to XYZ coordinates, with the origin at the center of the image
    """
    import numpy as np
    coordinates = np.array(coordinates)
    center = np.array(center)
    transformed_coordinates = np.array([(center[2]+p[2], center[1]-p[1], center[0]+p[0]) for p in coordinates])
    
    return np.clip(transformed_coordinates, -shape[0]//2, shape[0]//2-1)
    
    
def transform_coordinates_using_rotation_matrix(coordinates, rotation_matrix):
    """
    Transform a set of XYZ coordinates using a rotation matrix
    """
    import numpy as np
    coordinates = np.array(coordinates)
    coordinates = np.dot(coordinates, rotation_matrix)
    return coordinates 

def trilinear_interpolation(FT, rotation_matrix):
    """
    Perform trilinear interpolation on a 3D FT image
    """
    import numpy as np
    import numpy.linalg as la
    from scipy.interpolate import interpn
    
    rotation_matrix = np.transpose(rotation_matrix)
    ## FT = Fourrier transform of the map

    # Get the shape of the FT image
    nz, ny, nx = FT.shape

    depth, row, column = np.mgrid[:nz, :ny, :nx]
    numpy_coordinates = np.column_stack((depth.ravel(), row.ravel(), column.ravel()))  ## numpy_coordinates is ZYX format
    values_at_coordinates = FT.ravel()

    # Limit calculation to one half of a spherical mask within the FT

    # Get the spherical mask
    radius_pixels = nx//2
    spherical_mask = get_spherical_mask(FT.shape, radius_pixels).astype(bool)

    # Make False the spherical mask values corresponding to negative X axis
    hemisphere_mask = spherical_mask.copy()
    x_min = 0
    x_max = nx//2
    hemisphere_mask[:,:,x_min:x_max] = False

    # Coordinates of hemisphere mask where True
    hemisphere_coordinates = np.asarray(np.where(hemisphere_mask)).T
    

    # Convert hemisphere coordinates to XYZ coordinates
    center = np.array([nx//2, ny//2, nz//2])
    xyz_hemisphere_coordinates = transform_origin_to_center_of_image(hemisphere_coordinates, center, FT.shape)
    print("Transformed hemisphere coordinates")
    
    h, k, l = xyz_hemisphere_coordinates.T.astype(int)
    hT, kT, lT = -1*h, -1*k, -1*l
    xyz_hemisphere_coordinates_conjugate = -1 * xyz_hemisphere_coordinates
    conjugate_coordinates = transform_coordinates_to_numpy_coordinates(xyz_hemisphere_coordinates_conjugate, center, FT.shape).astype(int)
    z_conj, y_conj, x_conj = conjugate_coordinates.T
    conjugate_mask = np.zeros(FT.shape)
    conjugate_mask[z_conj, y_conj, x_conj] = 1
    conjugate_mask = conjugate_mask.astype(bool)
    

    # Transform XYZ hemisphere coordinates using rotation matrix
    xyz_hemisphere_coordinates_transform = transform_coordinates_using_rotation_matrix(xyz_hemisphere_coordinates, rotation_matrix)
    print("Transformed XYZ hemisphere coordinates")

    # Interpolation points in numpy coordinates
    #interpolation_points = transform_coordinates_to_numpy_coordinates(xyz_hemisphere_coordinates_transform, center, FT.shape)
   # print("Transformed interpolation points")#

    # Interpolate
    zlims = np.linspace(int(-nz//2), int(nz//2)-1, nz)
    ylims = np.linspace(int(-ny//2), int(ny//2)-1, ny)
    xlims = np.linspace(int(-nx//2), int(nx//2)-1, nx)
    
    interpolated_values = interpn((zlims, ylims, xlims), FT, xyz_hemisphere_coordinates_transform, method='linear')
    print("Interpolated values")

    # Apply interpolated values to a new FT image
    interpolated_FT = np.zeros(FT.shape,dtype="complex_")
    interpolated_FT[hemisphere_mask] = interpolated_values

    interpolated_FT[lT,kT,hT] = np.conjugate(interpolated_FT[l,k,h])
    print("is nan? ",np.isnan(interpolated_FT).any())

    return interpolated_FT
    



    



    





































# # Spherical mask on the FT image
#     spherical_indices = get_spherical_mask(shape=FT.shape, radius_pixels=nx//2).astype(bool)

#     # Get the coordinates of all non zero elements in the spherical mask
    
#     numpy_coordinates = np.asarray(np.where(spherical_indices)).T
#     values_in_spherical_mask = FT[spherical_indices] ## numpy coordinates are in the form of (z,y,x)

#     ## Convert numpy coordinates to XYZ coordinates
#     # Get the center of the image
#     center = np.array([nz//2, ny//2, nx//2])
#     # Transform the numpy coordinates to XYZ coordinates
#     LKH_coordinates_full = transform_coordinates_to_center_of_image(numpy_coordinates, center)
   
#     # Make values of mask corresponding to negative X frequencies zero
#     spherical_indices_halfmask = spherical_indices.copy()
#     x_min = 0
#     x_max = nx//2
#     spherical_indices_halfmask[:,:,x_min:x_max] = False

#     # Obtain coordinates of all non-zero pixels in the spherical mask
#     numpy_coordinates_halfmask = np.asarray(np.where(spherical_indices_halfmask)).T
           
#     # Transform the coordinates of the spherical mask to XYZ coordinates
#     LKH_coordinates_halfmask = transform_coordinates_to_center_of_image(numpy_coordinates_halfmask, center)

#     # Find the symmetrically corresponding coordinates in the FT image using the rotation matrix
#     LKH_symmetrical_coordinates = transform_coordinates_using_rotation_matrix(LKH_coordinates_halfmask, np.transpose(rotation_matrix))

#     # Obtain structure factors at the symmetrically corresponding coordinates using trilinear interpolation
#     values_in_halfmask = interpn(points=LKH_coordinates_full, values=values_in_spherical_mask, xi=LKH_symmetrical_coordinates, bounds_error=False, fill_value=0)

#     # Combine the values of the spherical mask with the values of the structure factors
#     new_FT = np.zeros(FT.shape)
#     new_FT[spherical_indices_halfmask] = values_in_halfmask

#     invert_spherical_indices = np.logical_not(spherical_indices_halfmask)
#     new_FT[invert_spherical_indices] = np.conj(np.flip(new_FT[spherical_indices_halfmask], axis=2))

#     return new_FT


    






    



    