#!/usr/bin/env python 
"""SWMatch - a Python class
"""

__author__ = "Longbin Chen(longbin@yahoo-inc.com, lbchen@cs.ucsb.edu)"
__version__ = "$Revision: 2.0 $"
__history__="""

    """
import sys
from string import *
import os
from math import *
import glob
from scipy import * 
import math
import getopt
import fileinput
import string   
import pickle
from opencv import cv
from opencv import highgui
from CntPoint import *
from CntAngle import *
from CntDAS import *
from SW import *


OP_MATCH = 1
OP_SUBST = 2
OP_DELET = 3 
OP_INSRT = 4 
NO_OP =  -1
NO_OP_BOTH = -1
NO_OP_A =   -2
NO_OP_B = -3 
  
 

def feature_filter(featurenames):
    	res = []
    	for f in featurenames:
    	    if (not(f in ['index','x','y'])):
    	    	res.append(f)
    	return res
        
def getdata(keypoints):
        res = keypoints[0].getvalues().keys()
        res.sort()
        feature = []
        featurenames = feature_filter(res)
        for pt in keypoints:
            val = []
            var = pt.getvalues()
            for i in range(len(featurenames)):
                val.append(var[featurenames[i]])
            feature.append(val)
        return feature


def usage():
    print "\nMatching two images using SW.py"
    print "============================================="
    print "usage: %s  [-o outputfile]  [-d] [-n pointnum] [-e] imagefile1 imagefile2"%sys.argv[0]
    print "Options:"
    print "    -e: even sampling"
    print "    -d: display result"
    print "    -o: output file name"
    print "    -h: help"
    print "    -n: feature point number"
    print "==========Dependency==================================="
    print "CntPoint.py, CntDAS.py CntAngle.py SW.py"
    print "======================================================="
    print "longbin chen, longbinc@yahoo.com"    

def main():

    ct1 = CntPoint()
    ct2 = CntPoint()
    agl = CntAngle()
    das = CntDAS()
    das.bDraw = 0


    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:dn:e", ["help", "output=", "draw", "num=", "even"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    output = None
    bDraw = 0
    npoint = 100
     
    for o, a in opts:
        if o == "-v":
            ct.verbose = 1
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--output"):
            output = a
        if o in ("-d", "--draw"):
            bDraw = 1
        if o in ("-n", "--num"):
            npoint = string.atoi(a)
        if o in ("-e", "--even"):
            ct.bEven = 1
    if (len(args)) != 2:
        usage()
        sys.exit(2)
    
    ct1.GetContour(args[0], npoint)
    
    agl.ExtractFeature(ct1.GetKeyPoints(), ct1.drawimg)
    das.ExtractFeature(ct1.GetKeyPoints(), ct1.drawimg)
    
    ct2.GetContour(args[1], npoint)
    agl.ExtractFeature(ct2.GetKeyPoints(), ct2.drawimg)
    das.ExtractFeature(ct2.GetKeyPoints(), ct2.drawimg)

    
    seq1 = getdata(ct1.GetKeyPoints())
    seq2 = getdata(ct2.GetKeyPoints())
    matcher = SmithWaterman()
    cost,align,X,Y = matcher.Align(seq1, seq2)
    myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
    if (bDraw):
        ct1.DrawKeyPoints()
        kpoints1 = ct1.GetKeyPoints()
        ct2.DrawKeyPoints()
        kpoints2 = ct2.GetKeyPoints()
        ptcount = 0
        for i in range(len(X)):
            xi = X[i]
            yi = Y[i]
            if (xi == -1):
                cv.cvPutText(ct2.drawimg, 'O', cv.cvPoint(int(kpoints2[yi].x), int(kpoints2[yi].y)), myfont, cv.cvScalar(255, 0, 0,0))
            if (yi == -1):
                cv.cvPutText(ct1.drawimg, 'O', cv.cvPoint(int(kpoints1[xi].x), int(kpoints1[xi].y)), myfont, cv.cvScalar(255, 0, 0,0))
            if (xi != -1 and yi != -1):
                ptcount  += 1
                cv.cvPutText(ct1.drawimg, str(ptcount), cv.cvPoint(int(kpoints1[xi].x), int(kpoints1[xi].y)), myfont, cv.cvScalar(255, 255, 0,0))
                cv.cvPutText(ct2.drawimg, str(ptcount), cv.cvPoint(int(kpoints2[yi].x), int(kpoints2[yi].y)), myfont, cv.cvScalar(255, 255, 0,0))
            
        highgui.cvNamedWindow ("contour1", 1)
        highgui.cvNamedWindow ("contour2", 1)
        highgui.cvShowImage ("contour1", ct1.drawimg)
        highgui.cvShowImage ("contour2", ct2.drawimg)
        highgui.cvWaitKey (0)       

     	
if __name__ == '__main__':
    # load the image gived on the command line
    main()
