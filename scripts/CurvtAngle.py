#! /usr/bin/env python
 

# import the necessary things for OpenCV

import sys
import math

import getopt
import string
import fileinput
from opencv import cv
from opencv import highgui
from CntPoint import *
 
  

class CurveAngle:
 
        
        
def usage():
    print "==========Usage========================================"
    print "usage: %s [-h] [-o outputfile]  [-d]  imagefile"%sys.argv[0]
    print "Options:"
    print "    -o: output file name"
    print "    -h: help"
    print "==========Dependency==================================="
    print "CntPoint.py CntAngle.py"
    print "======================================================="
    print "longbin chen, longbinc@yahoo.com"
    
def main():

    ct = CurvePoint()
    agl = CntAngle()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:", ["help", "output="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    output = None
     
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--output"):
            output = a

    if (len(args)) != 1:
        usage()
        sys.exit(2)

    for c in ct.allselected:
        agl.ExtractFeature(c, 0)

    if (output):
        ct.Save(output)

    
if __name__ == '__main__':
    # load the image gived on the command line
    main()
