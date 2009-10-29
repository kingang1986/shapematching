#!/usr/bin/env python 
 
import sys, math, os, optparse, string, fileinput
from MSS import *
from exDAS import *
from exAngle import *
from exSC import *

 
   
def main():

    usage = "%prog [options] <msshape> <outputmss>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)

    das = exDAS()
    agl = exAngle()
    sc = exSC()
    mss = MSS()
    mss.load(args[1])
    for seq in mss.seqs:
        das.ExtractFeature(seq.points, None)
        agl.ExtractFeature(seq.points, None)
    sc.ExtractFeature(mss)
    mss.save(args[2])
    

if __name__ == '__main__':
    # load the image gived on the command line
    main()
 
