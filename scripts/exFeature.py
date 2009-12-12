#!/usr/bin/env python 
 
import sys, math, os, optparse, string, fileinput
from MSS import *
from exDAS import *
from exAngle import *
from exSC import *
from exOffset import *


 
   
def main():

    usage = "%prog [options] <msshape> <outputmss>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-a', '--angle', dest = 'angle', default = False, action = 'store_true', help = 'generate angle feature')
    oparser.add_option('-d', '--das', dest = 'das', default = False, action = 'store_true', help = 'generate das feature')
    oparser.add_option('-s', '--sc', dest = 'sc', default = False, action = 'store_true', help = 'generate shape context feature')
    oparser.add_option('-o', '--offset', dest = 'offset', default = False, action = 'store_true', help = 'generate offset feature')


    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    if (not options.angle and not options.das and not options.sc and not options.offset):
        options.angle = True
        options.das   = True
        options.sc    = True
        options.offset= True

    das = exDAS()
    agl = exAngle()
    sc = exSC()
    mss = MSS()
    mss.load(args[1])
    for seq in mss.seqs:
        if (options.das):
            das.ExtractFeature(seq.points, None)
        if (options.angle):
            agl.ExtractFeature(seq.points, None)
    if (options.sc):
        keypoints = mss.getAllPoints()
        sc.ExtractFeature(keypoints, None)
    if (options.offset):
        oft = exOffset()
        oft.ExtractFeature(keypoints, None)
    mss.save(args[2])
    

if __name__ == '__main__':
    # load the image gived on the command line
    main()
 
