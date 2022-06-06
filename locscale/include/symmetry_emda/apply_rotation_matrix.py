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
    Transform a set of 3D numpy indices to hkl coordinates, with the origin at the center of the image
    """
    import numpy as np
    coordinates = np.array(coordinates)
    center = np.array(center)
    transformed_coordinates = coordinates - center  ## still in ZYX
    # convert to hkl
    hkl_transformed_coordinates = np.flip(transformed_coordinates, axis=-1) ## now XYZ > HKL
    hkl_transformed_coordinates[:,1] *= -1
    #print((hkl_transformed_coordinates[:,2]==transformed_coordinates[:,0]).all())
    #print((hkl_transformed_coordinates[:,0]==transformed_coordinates[:,2]).all())
    # Invert y axis
    
    return hkl_transformed_coordinates
    #return np.clip(invert_y_axis, 0, shape[0]//2-1)


def transform_coordinates_to_numpy_coordinates(coordinates, center, shape):
    """
    Transform a set of 3D numpy indices to hkl coordinates, with the origin at the center of the image
    """
    import numpy as np
    coordinates = np.array(coordinates)
    center = np.array(center)
    transformed_coordinates = np.array([(center[2]+p[2], center[1]-p[1], center[0]+p[0]) for p in coordinates])
    
    return np.clip(transformed_coordinates, -shape[0]//2, shape[0]//2-1)
    
    
def transform_coordinates_using_rotation_matrix(coordinates, rotation_matrix):
    """
    Transform a set of hkl coordinates using a rotation matrix
    """
    import numpy as np
    coordinates = np.array(coordinates)
    coordinates_T = np.matmul(rotation_matrix, np.transpose(coordinates))
    coordinates_T = np.transpose(coordinates_T)
        
    
    return coordinates_T 

def get_depth_row_col_ix_point(x,center):
    nz, ny, nx = center
    
    col_ix = int(x[0]+nx//2)
    row_ix = int(ny//2 - x[1])
    depth_ix = int(x[2]+nz//2)
    
    return depth_ix, row_ix, col_ix

def get_depth_row_col_ix(x,center):
    nz, ny, nx = center
    
    col_ix = x[:,0]+nx//2
    row_ix = ny//2 - x[:,1]
    depth_ix = x[:,2]+nz//2
    
    return depth_ix.astype(int), row_ix.astype(int), col_ix.astype(int)
    
def trilinear_interpolation(FT, rotation_matrix):
    import numpy as np
    from tqdm import tqdm
    from scipy.interpolate import RegularGridInterpolator
    #from scipy.iinterpolate import interpn
    
    
    nz, ny, nx = FT.shape
    radius_pixels = nx//2
    spherical_mask = get_spherical_mask(FT.shape, radius_pixels).astype(bool)
    interpolate_FT = np.zeros(FT.shape, dtype="complex_")
    
    zlims = np.linspace(int(-nz//2), int(nz//2)-1, nz)
    ylims = np.linspace(int(-ny//2), int(ny//2)-1, ny)
    xlims = np.linspace(int(-nx//2), int(nx//2)-1, nx)
    
    interpolator = RegularGridInterpolator(points=(zlims, ylims, xlims), \
                                                  values=FT, \
                                                  method='linear')
    
    rotation_matrix_transpose = np.transpose(rotation_matrix)
    center = (nz//2, ny//2, nx//2)
    pbar = tqdm(total=nx//2*ny*nz, desc="Symmetrising: ")
    slist = []
    cdlist = []
    conjlist = []
    compare_interpolation_methods = []
    for depth in range(nz):
        for row in range(ny):
            for col in range(nx//2):
                pbar.update(1)
                if spherical_mask[depth,row,col]: 
                        
                    h = col-nx//2 ## center
                    k = -1 * (row-ny//2) ## flip Y for correct axis
                    l = depth - nz//2 
                    slist.append(np.array([h,k,l]))
                    cdlist.append(np.array([depth, row, col]))
                    
                    s = np.array([h,k,l])
                    x = np.matmul(rotation_matrix_transpose, s)
                    x0 = np.floor(x)
                    x1 = x0+1
                    xd = x-x0
                    xdr = 1 - xd
                    assert (xd<1).all(), "xd: {}".format(xd)
                    
                    d0, r0, c0 = get_depth_row_col_ix_point(x0, center)
                    d1, r1, c1 = get_depth_row_col_ix_point(x1, center)
                    
                                        
                    c000 = FT[d0,r0,c0]
                    c001 = FT[d0,r0,c1]
                    c010 = FT[d0,r1,c0]
                    c011 = FT[d0,r1,c1]
                    c100 = FT[d1,r0,c0]
                    c101 = FT[d1,r0,c1]
                    c110 = FT[d1,r1,c0]
                    c111 = FT[d1,r1,c1]
                    
                    c00 = c000*xdr[0] + c100*xd[0]
                    c01 = c001*xdr[0] + c101*xd[0]
                    c10 = c010*xdr[0] + c110*xd[0]
                    c11 = c011*xdr[0] + c111*xd[0]
                    
                    c0 = c00*xdr[1] + c10*xd[1]
                    c1 = c01*xdr[1] + c11*xd[1]
                    
                    c = c0*xdr[2] + c1*xd[2]
                    
                    c_scipy = interpolator(np.flip(x))
                    #compare_interpolation_methods.append(c==c_scipy)
                    if c != c_scipy:
                        print("not equal")
                        print(c)
                        print(c_scipy)
                        print(h,k,l)
                    if c == c_scipy:
                        print("equal")
                        print("here > ",h,k,l)
                    
                    dC, rC, cC = get_depth_row_col_ix_point(-1*s, center)
                    conjlist.append(np.array([dC,rC,cC]))
                    interpolate_FT[depth, row, col] = c
                    interpolate_FT[dC, rC, cC] = np.conjugate(c)
                    
                                        
                else:
                    continue
    
    slist = np.array(slist)
    print("smin" , slist.min())
    print("smax" , slist.max())
    cdlist = np.array(cdlist)
    conjlist = np.array(conjlist)
    compare_interpolation_methods = np.array(compare_interpolation_methods)
    
    test_hkl_transform = transform_origin_to_center_of_image(cdlist, center, FT.shape)
    
    dtest, rtest, ctest = get_depth_row_col_ix(-1*slist, center)
    
    
    print("coordinate check: ",(test_hkl_transform==slist).all())
    
    print("conjugate coordinate check D: ", (dtest==conjlist[:,0]).all())
    print("conjugate coordinate check R: ", (rtest==conjlist[:,1]).all())
    print("conjugate coordinate check C: ", (ctest==conjlist[:,2]).all())
    
    print("interpolation alg check: ", compare_interpolation_methods.all())
    print("is nan? ",np.isnan(interpolate_FT).any())
        
    x_slice_FT = abs(FT[:,:,140])
    x_slice_interpolated_FT = abs(interpolate_FT[:,:,140])
    
    y_slice_FT = abs(FT[:,140,:])
    y_slice_interpolated_FT = abs(interpolate_FT[:,140,:])
    
    z_slice_FT = abs(FT[140,:,:])
    z_slice_interpolated_FT = abs(interpolate_FT[140,:,:])
    
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
    
    print("roptation: \n", rotation_matrix.round(2))
    
    rscc_slice_x = rsc(x_slice_FT, x_slice_interpolated_FT)
    rscc_slice_y = rsc(y_slice_FT, y_slice_interpolated_FT)
    rscc_slice_z = rsc(z_slice_FT, z_slice_interpolated_FT)
    
    print("rscc x slice:",rscc_slice_x)
    print("rscc y slice:",rscc_slice_y)
    print("rscc z slice:",rscc_slice_z)
    
    return interpolate_FT
                

def trilinear_interpolation2(FT, rotation_matrix):

    #Perform trilinear interpolation on a 3D FT image

    import numpy as np
    #import numpy.linalg as la
    from scipy.interpolate import interpn
    from scipy.interpolate import RegularGridInterpolator
    
    #rotation_matrix = np.transpose(rotation_matrix)
    ## FT = Fourrier transform of the map

    # Get the shape of the FT image
    nz, ny, nx = FT.shape

    #depth, row, column = np.mgrid[:nz, :ny, :nx]
    #numpy_coordinates = np.column_stack((depth.ravel(), row.ravel(), column.ravel()))  ## numpy_coordinates is ZYX format
    #values_at_coordinates = FT.ravel()

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
    hemisphere_coordinates = np.asarray(np.where(hemisphere_mask)).T  # ZYX
    

    # Convert hemisphere coordinates to hkl coordinates
    center = np.array([nx//2, ny//2, nz//2])
    hkl_hemisphere_coordinates = transform_origin_to_center_of_image(hemisphere_coordinates, center, FT.shape)
    print("Transformed hemisphere coordinates")
    
    
    
#    hkl_hemisphere_coordinates_conjugate = -1 * hkl_hemisphere_coordinates
#    conjugate_coordinates = transform_coordinates_to_numpy_coordinates(hkl_hemisphere_coordinates_conjugate, center, FT.shape).astype(int)
#    z_conj, y_conj, x_conj = conjugate_coordinates.T
#    conjugate_mask = np.zeros(FT.shape)
#    conjugate_mask[z_conj, y_conj, x_conj] = 1
#    conjugate_mask = conjugate_mask.astype(bool)
    

    # Transform hkl hemisphere coordinates using rotation matrix
    hkl_hemisphere_coordinates_transform = transform_coordinates_using_rotation_matrix(hkl_hemisphere_coordinates, rotation_matrix)
    hkl_hemisphere_coordinates_transform = np.clip(hkl_hemisphere_coordinates_transform, int(-nz//2), int(nz//2)-1)
    print("Transformed hkl hemisphere coordinates")
    #print(hkl_hemisphere_coordinates_transform[:,0].min(), hkl_hemisphere_coordinates_transform[:,0].max())
    #print(hkl_hemisphere_coordinates_transform[:,1].min(), hkl_hemisphere_coordinates_transform[:,1].max())
    #print(hkl_hemisphere_coordinates_transform[:,2].min(), hkl_hemisphere_coordinates_transform[:,2].max())
    
    

    # Interpolation points in numpy coordinates
    #interpolation_points = transform_coordinates_to_numpy_coordinates(hkl_hemisphere_coordinates_transform, center, FT.shape)
   # print("Transformed interpolation points")#

    # Interpolate
    zlims = np.linspace(int(-nz//2), int(nz//2)-1, nz)
    ylims = np.linspace(int(-ny//2), int(ny//2)-1, ny)
    xlims = np.linspace(int(-nx//2), int(nx//2)-1, nx)
    
    interpolator = RegularGridInterpolator(points=(zlims, ylims, xlims), \
                                                  values=FT, \
                                                  method='linear')
    interpolated_values = interpolator(hkl_hemisphere_coordinates_transform)
    print("Interpolated values")

    # Apply interpolated values to a new FT image
    interpolated_FT = np.zeros(FT.shape,dtype="complex_")
    interpolated_FT[hemisphere_mask] = interpolated_values
    depth, row, col = get_depth_row_col_ix(hkl_hemisphere_coordinates, center)
    depth_T, row_T, col_T = get_depth_row_col_ix(-1*hkl_hemisphere_coordinates, center)
    
    interpolated_FT[depth_T,row_T,col_T] = np.conjugate(interpolated_FT[depth,row,col])
    print("is nan? ",np.isnan(interpolated_FT).any())
        
    x_slice_FT = abs(FT[:,:,140])
    x_slice_interpolated_FT = abs(interpolated_FT[:,:,140])
    
    y_slice_FT = abs(FT[:,140,:])
    y_slice_interpolated_FT = abs(interpolated_FT[:,140,:])
    
    z_slice_FT = abs(FT[140,:,:])
    z_slice_interpolated_FT = abs(interpolated_FT[140,:,:])
    
    from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rsc
    
    print("roptation: \n", rotation_matrix.round(2))
    
    rscc_slice_x = rsc(x_slice_FT, x_slice_interpolated_FT)
    rscc_slice_y = rsc(y_slice_FT, y_slice_interpolated_FT)
    rscc_slice_z = rsc(z_slice_FT, z_slice_interpolated_FT)
    
    print("rscc x slice:",rscc_slice_x)
    print("rscc y slice:",rscc_slice_y)
    print("rscc z slice:",rscc_slice_z)
    return interpolated_FT
   


    



    





































# # Spherical mask on the FT image
#     spherical_indices = get_spherical_mask(shape=FT.shape, radius_pixels=nx//2).astype(bool)

#     # Get the coordinates of all non zero elements in the spherical mask
    
#     numpy_coordinates = np.asarray(np.where(spherical_indices)).T
#     values_in_spherical_mask = FT[spherical_indices] ## numpy coordinates are in the form of (z,y,x)

#     ## Convert numpy coordinates to hkl coordinates
#     # Get the center of the image
#     center = np.array([nz//2, ny//2, nx//2])
#     # Transform the numpy coordinates to hkl coordinates
#     LKH_coordinates_full = transform_coordinates_to_center_of_image(numpy_coordinates, center)
   
#     # Make values of mask corresponding to negative X frequencies zero
#     spherical_indices_halfmask = spherical_indices.copy()
#     x_min = 0
#     x_max = nx//2
#     spherical_indices_halfmask[:,:,x_min:x_max] = False

#     # Obtain coordinates of all non-zero pixels in the spherical mask
#     numpy_coordinates_halfmask = np.asarray(np.where(spherical_indices_halfmask)).T
           
#     # Transform the coordinates of the spherical mask to hkl coordinates
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


    






    



    
