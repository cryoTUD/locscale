import numpy as np

try:
    import pyfftw
    pyfftw_flag = True
except ImportError:
    pyfftw_flag = False


def plan_fft(arr,keep_shape=False,new_inparray=False):
    input_dtype = str(arr.dtype)
    if not keep_shape:
        output_dtype = 'complex64'
        if not input_dtype in ['float32','float64','longdouble']: 
            input_dtype = 'float32'
        elif input_dtype == 'float64':
            output_dtype = 'complex128'
        elif input_dtype == 'longdouble':
            output_dtype = 'clongdouble'
        #for r2c transforms:
        output_array_shape = \
                arr.shape[:len(arr.shape)-1]+ \
                                    (arr.shape[-1]//2 + 1,)
    else:
        output_dtype = 'complex64'
        output_array_shape = arr.shape
    
    fftoutput = pyfftw.n_byte_align_empty(output_array_shape, 
                                    n=16, dtype=output_dtype)
    #check if array is byte aligned
    #TODO: can we read the map file as byte aligned? 
    if new_inparray or not pyfftw.is_byte_aligned(arr):
        inputarray = pyfftw.empty_aligned(arr.shape,
                                          n=16,dtype='float32')
        fft = pyfftw.FFTW(inputarray,fftoutput,
                          direction='FFTW_FORWARD',axes=(0,1,2),
                          flags=['FFTW_ESTIMATE'])
    elif pyfftw.is_byte_aligned(arr):
        fft = pyfftw.FFTW(arr,fftoutput,
                          direction='FFTW_FORWARD',axes=(0,1,2),
                          flags=['FFTW_ESTIMATE'])
        inputarray = arr
        
    return fft, fftoutput, inputarray

def calculate_fft(arr,keep_shape=False):
    
    if pyfftw_flag:
        fft, fftoutput, inputarray = plan_fft(arr, keep_shape=keep_shape)
        inputarray[:,:,:] = arr
        fft()
    else:
        #TODO: warning raises error in tasks
        #warnings.warn("PyFFTw not found!, using numpy fft")
        if not keep_shape:
            fftoutput = np.fft.rfftn(arr) 
        else:
            fftoutput = np.fft.fftn(arr)
    return fftoutput

def plan_ifft(arr,output_shape=None,output_array_dtype=None,
              new_inparray=False):
    input_dtype = str(arr.dtype)
#         #for c2r transforms:
#             if output_shape is None: output_shape = \
#                                     arr.shape[:len(arr.shape)-1]+\
#                                     ((arr.shape[-1] - 1)*2,)
    if output_array_dtype is None: output_array_dtype = 'float32'
    if output_shape is None: 
        output_shape = arr.shape[:len(arr.shape)-1]+\
                        ((arr.shape[-1] - 1)*2,)
        if not input_dtype in ['complex64','complex128','clongdouble']: 
            input_dtype = 'complex64'
        elif input_dtype == 'complex128':
            output_array_dtype = 'float64'
        elif input_dtype == 'clongdouble':
            output_array_dtype = 'longdouble'
    elif output_shape[-1]//2 + 1 == arr.shape[-1]:
        if not input_dtype in ['complex64','complex128','clongdouble']: 
            input_dtype = 'complex64'
        elif input_dtype == 'complex128':
            output_array_dtype = 'float64'
        elif input_dtype == 'clongdouble':
            output_array_dtype = 'longdouble'
    else:                           
        output_shape = arr.shape
        output_array_dtype = 'complex64'
        
    output_array = pyfftw.empty_aligned(output_shape, 
                                        n=16, dtype=output_array_dtype)
    #check if array is byte aligned
    if new_inparray or not pyfftw.is_byte_aligned(arr):
        inputarray = pyfftw.n_byte_align_empty(arr.shape,
                                               n=16,
                                               dtype=input_dtype)
        ifft = pyfftw.FFTW(inputarray,output_array,
                           direction='FFTW_BACKWARD',\
                       axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
        inputarray[:,:,:] = arr
    else:
        ifft = pyfftw.FFTW(arr,output_array,
                           direction='FFTW_BACKWARD',\
                           axes=(0,1,2),flags=['FFTW_ESTIMATE'])#planning_timelimit=0.5)
        inputarray = arr

    return ifft, output_array, inputarray

def calculate_ifft(arr,output_shape=None,inplace=False):
    """
    Calculate inverse fourier transform
    """
    if pyfftw_flag:
        ifft, output_array, inputarray = plan_ifft(arr,output_shape=output_shape)
        #r2c fft
        ifft()
    else:
        #TODO: warnings raises error in tasks
        #warnings.warn("PyFFTw not found!, using numpy fft")
        if output_shape is None or output_shape[-1]//2 + 1 == arr.shape[-1]:
            output_array = np.fft.irfftn(arr)
        else:
            output_array = np.real(np.fft.ifftn(arr))
        #np.fft.fftpack._fft_cache.clear()    
    del arr
    return output_array.real.astype(np.float32,copy=False)

