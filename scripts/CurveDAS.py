#!/usr/bin/env python
 

# import the necessary things for OpenCV

import sys, math, optparse, fileinput
import math

from opencv import cv
from opencv import highgui
from CntPoint import *
from CurvePoint import *
from CntDAS import * 



def main():

    usage = "%prog [options] <msshape> <outputmss>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    #oparser.add_option('-c', '--shapecontext', action="store_false",dest = 'shapecontext', default = True,  help = 'extract shape context  features(default true)')
    #oparser.add_option('-w', '--das', action="store_false",dest = 'das', default = True,  help = 'extract DAS features(default true)')
    
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)

    das = CntDAS()
    mss = MSS()
    mss.load(args[1])
    for seq in mss.seqs:
        das.ExtractFeature(seq.points, None)
    mss.save(args[2])
    
if __name__ == '__main__':
    # load the image gived on the command line
    main()
