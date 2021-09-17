from EMAN2 import *
import os

cmdl_parser = argparse.ArgumentParser(
                prog="set_ori_zero.py",
                description='*** Make MRC map origin zero-based ***',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

required_args = cmdl_parser.add_argument_group("Input")
optional_args = cmdl_parser.add_argument_group("Options")

required_args.add_argument('-m', '--map', metavar="map.ccp4",  type=argparse.FileType('r'), required=True,
                         help='Input filename map')
optional_args.add_argument('-h', '--help', action="help",
                         help="show this help message and exit")

args = cmdl_parser.parse_args()


def set_zero_origin(input_map):
    """
    Make map origin zero-based
    """
    map = EMData()
    map.read_image(os.path.abspath(args.map.name),0,True)

    print "Map origin was:\t\t",map['MRC.nxstart'],map['MRC.nystart'],map['MRC.nzstart']
    map['MRC.nxstart'] = 0
    map['MRC.nystart'] = 0
    map['MRC.nzstart'] = 0
    print "Map origin set to:\t",map['MRC.nxstart'],map['MRC.nystart'],map['MRC.nzstart']

    map.write_image(os.path.abspath(args.map.name),0,EMUtil.ImageType.IMAGE_MRC, True)


if (__name__ == "__main__") :
   set_zero_origin(args.map)

