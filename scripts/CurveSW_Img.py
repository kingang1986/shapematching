#!/home/longbin/install/apps/Python-2.5_gcc/bin/python
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
from CntSC import *
from SW import *
from CurvePoint import *
from CurveAngle import *
from myroutine import *

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
def putoriginal(fname, img):
    ori_img = highgui.cvLoadImage (fname)
    ori_img_thumb = cv.cvCreateImage(cv.cvSize(ori_img.width/4, ori_img.height/4), 8,3)
    cv.cvResize(ori_img, ori_img_thumb)
    for x in range(ori_img_thumb.height):
        for y in range(ori_img_thumb.width):
            cv.cvSet2D(img, x, y, cv.cvGet2D(ori_img_thumb, x, y))
    return 
    
     


def main():

    ct1 = CurvePoint()
    ct2 = CurvePoint()
    sc = CntSC()
    ang = CntAngle()
     
    
    usage = "%prog [options] <imgfile1> <imgfile2>"
    version = "%prog 0.2\nLongbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--display', action="store_true", dest = 'display', default = False, help = 'display the image')
    oparser.add_option('-n', '--number', dest = 'num',  type="int", default = 200 , help = 'the number of feature points')
    oparser.add_option('-s', '--save', dest = 'save', default = None, help = 'save the img file')

    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)

    ct1.GetContour(args[1], options.num)
    allkeys = []
    for c in ct1.allselected:
        allkeys = allkeys + c
    sc.ExtractFeature(allkeys)
    ang.ExtractFeature(allkeys,0); 
    allkeys = []
    ct2.GetContour(args[2], options.num)
    for c in ct2.allselected:
        allkeys = allkeys + c
    sc.ExtractFeature(allkeys)
    ang.ExtractFeature(allkeys,0); 

    sumscore = []
    matcher = SmithWaterman()
    ct1.bDrawNumber = 0
    ct2.bDrawNumber = 0
    if (options.display):
        ct1.DrawKeyPoints()
        ct2.DrawKeyPoints()
    myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
    idx = -1
    putoriginal(args[1], ct1.drawimg)
    putoriginal(args[2], ct2.drawimg)
    cv.cvNot(ct1.drawimg, ct1.drawimg)
    cv.cvNot(ct2.drawimg, ct2.drawimg)
    for c1 in ct1.allselected:
        idx += 1
        cscore = -100000000
        cpt1 =   getdata(c1)
        bX = []
        bY = []
        bestcurve = None
        for c2 in ct2.allselected:
            cpt2 =   getdata(c2)
            cost,align,X,Y = matcher.Align(cpt1, cpt2)
            normalized_score = cost - log10(len(c2) + 1) * 1000
            print len(c1), len(c2),cost, normalized_score, cscore
            if (normalized_score > cscore):
                cscore = normalized_score
                bX = X[:]
                bY = Y[:]
                bestcurve = c2
        if (options.display):
            ptcount = 0
            for i in range(len(bX)):
                xi = bX[i]
                yi = bY[i]
                #if (xi == -1):
                    #cv.cvDrawCircle(ct2.drawimg, cv.cvPoint(int(bestcurve[yi].x), int(bestcurve[yi].y)),4, cv.cvScalar(255,0,0,0))
                    #cv.cvPutText(ct2.drawimg, 'O', cv.cvPoint(int(c2[yi].x), int(c2[yi].y)), myfont, cv.cvScalar(255, 0, 0,0))
                #if (yi == -1):
                    #cv.cvDrawCircle(ct1.drawimg, cv.cvPoint(int(c1[xi].x), int(c1[xi].y)),4, cv.cvScalar(255,0,0,0))
                    #cv.cvPutText(ct1.drawimg, 'O', cv.cvPoint(int(c1[xi].x), int(c1[xi].y)), myfont, cv.cvScalar(255, 0, 0,0))
                if (xi != -1 and yi != -1):
                    ptcount  += 1
                    cv.cvDrawCircle(ct1.drawimg, cv.cvPoint(int(c1[xi].x), int(c1[xi].y)),2, clrs[idx])
                    cv.cvPutText(ct1.drawimg, str(ptcount), cv.cvPoint(int(c1[xi].x), int(c1[xi].y)), myfont, clrs[idx])
                    cv.cvDrawCircle(ct2.drawimg, cv.cvPoint(int(bestcurve[yi].x), int(bestcurve[yi].y)),2, clrs[idx])
                    cv.cvPutText(ct2.drawimg, str(ptcount), cv.cvPoint(int(bestcurve[yi].x), int(bestcurve[yi].y)), myfont, clrs[idx])
        sumscore.append(cscore)
    print sumscore
    if (options.display):            
	    highgui.cvNamedWindow ("contour1", 1)
	    highgui.cvNamedWindow ("contour2", 1)
	    highgui.cvShowImage ("contour1", ct1.drawimg)
	    highgui.cvShowImage ("contour2", ct2.drawimg)
	    highgui.cvWaitKey (0)       
    if (options.save):
        mergeimg = mergeimage_83(ct1.drawimg, ct2.drawimg)
        highgui.cvSaveImage("_sw_result.bmp", mergeimg)
        #highgui.cvSaveImage("_sw_save1.bmp", ct1.drawimg)
        #highgui.cvSaveImage("_sw_save2.bmp", ct2.drawimg)
                    

        
if __name__ == '__main__':
    # load the image gived on the command line
    main()
