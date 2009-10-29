#!/home/longbin/install/apps/Python-2.5_gcc/bin/python 

# import the necessary things for OpenCV

import sys
import math

import optparse
import string
import fileinput
from opencv import cv
from opencv import highgui
from CurvePoint import *
from CntSC import *
    
def main():

    ct = CurvePoint()
    sc = CntSC()

    usage = "%prog [options] <pointfile>"
    version = "%prog 0.2\nLongbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 2:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)

    ct.LoadCont(args[1])
    allkeys = []
    for c in ct.allselected:
        allkeys = allkeys + c
    sc.ExtractFeature(allkeys)

    if (options.output):
        ct.Save(options.output)

    
if __name__ == '__main__':
    # load the image gived on the command line
    main()
