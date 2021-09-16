from EMAN2 import EMData, EMNumPy, Util, XYData
import numpy as np
from sparx import get_image, binarize, dilation, filt_gaussl, model_blank
import argparse, os, sys


progname = os.path.basename(sys.argv[0])
revision = filter(str.isdigit, "$Revision: 1 $")
datmod = "$Date: 2017-03-06 22:14:31 +0200 (Mo, 06 Mar 2017) $"      
author = 'authors: Arjen J. Jakobi' + '; ' + datmod [8:18]
version = progname + '  0.1' + '  (r' + revision + ';' + datmod [7:18] + ')'

simple_cmd = 'python automask_create.py -em emmap.mrc -p 1.0 -o mask.mrc'

cmdl_parser = argparse.ArgumentParser(
description='Example usage: \"{0}\"'.format(simple_cmd))

cmdl_parser.add_argument('-em', '--em_map', required=True, help='Input filename EM map')
cmdl_parser.add_argument('-p', '--apix', type=float, required=True, help='pixel size in Angstrom')
cmdl_parser.add_argument('-o', '--outfile', required=True, help='Output filename')


def set_origin_and_pixel_size(map, mask, pixel_size):
    orix, oriy, oriz = map['MRC.nxstart'], map['MRC.nystart'], map['MRC.nzstart']
    mask['MRC.nxstart'], mask['MRC.nystart'], mask['MRC.nzstart'] = orix, oriy, oriz
    mask.set_attr("apix_x", pixel_size)
    mask.set_attr("apix_y", pixel_size)
    mask.set_attr("apix_z", pixel_size)
    return mask

def build_structural_mask_from_volume(reference_vol, pixel_size, sigma_factor=1.0, width_falloff=0.03):
    mask = EMData()
    xsize, ysize, zsize = reference_vol.get_xsize(), reference_vol.get_ysize(), reference_vol.get_zsize()
    mask.set_size(xsize, ysize, zsize)
    mask.to_zero()
    if xsize == ysize and xsize == zsize and ysize == zsize:
        sphere_radius = mask.get_xsize() // 2
    else:
        sphere_radius = max(xsize, ysize, zsize) // 2
    mask.process_inplace("testimage.circlesphere", {"radius":sphere_radius})
    stat = Util.infomask(reference_vol, mask, True)
    binary_vol = binarize(reference_vol * mask, sigma_factor * stat[1])
    if width_falloff != 0:
            dilation_pixels = 7
            odd = 1
            if (dilation_pixels) % 2 != odd:
                    dilation_pixels += 1
            dilation_mask = model_blank(dilation_pixels, dilation_pixels, dilation_pixels, bckg = 1.0)
            structural_mask = dilation(binary_vol, dilation_mask)
    else:
            structural_mask = binary_vol
    return structural_mask

def create_auto_mask(args):
    map = get_image(args.em_map)
    automask = build_structural_mask_from_volume(map, args.apix)
    automask = set_origin_and_pixel_size(map, automask, args.apix)
    automask.write_image(args.outfile)

def main():
    args = cmdl_parser.parse_args()
    create_auto_mask(args)

if __name__ == '__main__':
    main()

