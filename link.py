#!/usr/bin/env python
import sys
import math

import string
import optparse
import fileinput
from CntPoint import *
from opencv import cv
from opencv import highgui
import random

_red = cv.cvScalar (0, 0, 255, 0)
_green = cv.cvScalar (0, 255, 0, 0)
_white = cv.cvScalar (255,255,255,0)
_black = cv.cvScalar (0,0,0,0)
color = [_red, _green, _white, cv.cvScalar(128,0, 128), cv.cvScalar(128,128, 128), cv.cvScalar(128,0, 0), cv.cvScalar(0,0, 128) ]

class Link:
        

    def __findedge(self, filename): #find the corners of images, and save all corner points in self.vKeyPoints
               
        tmpimg = highgui.cvLoadImage (filename)
        self.img = cv.cvCreateImage(cv.cvSize(tmpimg.width * 2, tmpimg.height * 2), 8, 3)
        cv.cvResize(tmpimg, self.img, cv.CV_INTER_LINEAR)
        #self.drawimg = cv.cvCloneImage(self.img)
        self.drawimg = cv.cvCreateImage(cv.cvGetSize(self.img), 8, 3)
        greyimg = cv.cvCreateImage(cv.cvSize(self.img.width, self.img.height), 8,1)
        cv.cvCvtColor (self.img, greyimg, cv.CV_BGR2GRAY)
        self.allcurve = []
        for i in range(80, 200, 40):
            bimg = cv.cvCloneImage(greyimg) 
            #cv.cvSmooth(bimg, bimg, cv.CV_MEDIAN, 9)
            cv.cvSmooth(bimg, bimg, cv.CV_MEDIAN, 9)
            cv.cvSmooth(bimg, bimg, cv.CV_BILATERAL, 9)
            cv.cvSmooth(bimg, bimg, cv.CV_BLUR, 9)
            cv.cvSmooth(bimg, bimg, cv.CV_BLUR, 9)
            cv.cvThreshold(greyimg, bimg, i, 255, cv.CV_THRESH_BINARY)
            self.__findcurve(bimg)
        
            

    def __findcurve(self, img):
        storage = cv.cvCreateMemStorage(0)
        nb_contours, cont = cv.cvFindContours (img,
            storage,
            cv.sizeof_CvContour,
            cv.CV_RETR_LIST,
            cv.CV_CHAIN_APPROX_NONE,
            cv.cvPoint (0,0))
        cidx = int(random.random() * len(color))
        #cv.cvDrawContours (self.drawimg, cont, _white, _white, 1, 1, cv.CV_AA, cv.cvPoint (0, 0))
        idx = 0
        for c in cont.hrange():
            PointArray = cv.cvCreateMat(1, c.total, cv.CV_32SC2)
            PointArray2D32f= cv.cvCreateMat( 1, c.total  , cv.CV_32FC2)
            cv.cvCvtSeqToArray(c, PointArray, cv.cvSlice(0, cv.CV_WHOLE_SEQ_END_INDEX))
            fpoints = []
            for i in range(c.total):
                kp = KeyPoint()
                kp.x = cv.cvGet2D(PointArray,0, i)[0]
                kp.y = cv.cvGet2D(PointArray,0, i)[1]
                kp.index = idx
                idx += 1
                fpoints.append(kp)
            self.allcurve.append(fpoints)
        self.curvelength = idx

    def __link(self):
        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        kkk = 0
        for curve in self.allcurve:
            showpt = []
            state = 0  
            currentPoint = None 
            cumulate = 0
            dcurve = curve + curve
            curlen = len(curve)
            ptcount = 0
            for c in dcurve:
                if (ptcount > curlen): break
                cumulate += 1
                for kk in range(len(self.points)):
                     k = self.points[kk]
                     if (abs(c.x - k.x) + abs(c.y - k.y) < 10):
                          #if (kk == 28):
                          #    print "hit curve %d state %d cum %d " % (len(curve), state, cumulate)
                          #    for t in curve:
                          #        cv.cvSet2D(self.drawimg, int(t.y), int(t.x), _green)
                          if (currentPoint != k or cumulate > 100):
                              state += 1
                              currentPoint = k
                              cumulate = 0
                if (state > 0):
                     showpt.append([c, state])
                     ptcount += 1
            if (state > 1):
                kkk += 1
                cnt = 0
                for s,t in showpt:
                    cnt += 1
                    if (t < state):
                         cv.cvSet2D(self.drawimg, int(s.y), int(s.x), color[kkk % 6 ])
                         #if (cnt / 40 * 40 == cnt):
                         #    cv.cvPutText(self.drawimg, str(kkk), cv.cvPoint(int(s.x), int(s.y)), myfont, _white)
        
    def DrawKeyPoints(self):
        if (not self.drawimg):
            self.drawimg = cv.cvCloneImage(self.img) 

        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        ic = 0
        for c in self.points:
            cv.cvPutText(self.drawimg, str(ic), cv.cvPoint(int(c.x), int(c.y)), myfont, cv.cvScalar(255, 255, 0,0))
            ic += 1
            cv.cvDrawCircle(self.drawimg, c, 4, cv.cvScalar(255,255,0,0))
        
    def SaveImage(self, filename):
        cv.cvNot(self.drawimg, self.drawimg)
        highgui.cvSaveImage(filename, self.drawimg)
 


    def LinkPoints(self, fname, ptfile, pointnum):
        self.nPointNum = pointnum
        self.__findedge(fname)
        self.__loadPoints(ptfile)
        self.__link()
        self.DrawKeyPoints()
    
    def __loadPoints(self, ptfile):
        self.points = []
        for line in fileinput.input(ptfile):
            dr = line.strip("\n").strip("\r").split(" ")
            ds = [d.strip(" ") for d in dr]
            dt = [d for d in dr if d != ""]
            x = float(dt[0]) * 2
            y = float(dt[1]) * 2
            pt = cv.cvPoint(int(x), int(y))
            self.points.append(pt)
              
    def __init__(self):
        self.output  = None
        self.bDraw  = 0
        self.bFilter = 1
        self.bDrawNumber = False
        self.drawimg = None
        
    
    
def main():
    ct = Link()
    usage = "%s [options]  <imgfile> <pointfile>" % (sys.argv[0])
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--display', action="store_true", dest = 'display', default = False, help = 'display the image')
    oparser.add_option('-n', '--number', dest = 'num', type='int',default = 200 , help = 'the number of feature points')
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')
    oparser.add_option('-s', '--save', dest = 'save', default = None, help = 'save the img file')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
   
    ct.LinkPoints(args[1], args[2], options.num)

    if (options.display):
        highgui.cvNamedWindow ("Corner1", 1)
        highgui.cvShowImage ("Corner1", ct.drawimg)
        highgui.cvWaitKey (0)   
    if (options.save):
        highgui.cvSaveImage(options.save, ct.drawimg)    

if __name__ == '__main__':
    # load the image gived on the command line
    main()
