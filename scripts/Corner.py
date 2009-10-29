#!/usr/bin/env python
import sys
import math

import string
import optparse
import fileinput
from CntPoint import *
from opencv import cv
from opencv import highgui

_red = cv.cvScalar (0, 0, 255, 0)
_green = cv.cvScalar (0, 255, 0, 0)
_white = cv.cvScalar (255,255,255,0)
_black = cv.cvScalar (0,0,0,0)

class Corner:

    def __FindCorner(self, filename): #find the corners of images, and save all corner points in self.vKeyPoints
        self.img = highgui.cvLoadImage (filename)
        greyimg = cv.cvCreateImage(cv.cvSize(self.img.width, self.img.height), 8,1)
        hsvimg = cv.cvCreateImage(cv.cvGetSize(self.img), 8, 3)
        cv.cvCvtColor(self.img, hsvimg, cv.CV_RGB2HSV)
        cv.cvCvtColor (hsvimg, greyimg, cv.CV_BGR2GRAY)
        
        eigImage = cv.cvCreateImage(cv.cvGetSize(greyimg), cv.IPL_DEPTH_32F, 1)
        tempImage = cv.cvCreateImage(cv.cvGetSize(greyimg), cv.IPL_DEPTH_32F, 1)
        self.points = cv.cvGoodFeaturesToTrack(greyimg, eigImage,tempImage, 2000, 0.01, 5, None, 3,0,0.01 )
        self.points2 = cv.cvFindCornerSubPix(greyimg, self.points,cv.cvSize(20, 20), 
                                             cv.cvSize(-1, -1), cv.cvTermCriteria(cv.CV_TERMCRIT_ITER |cv.CV_TERMCRIT_EPS, 20, 0.03))
        cv.cvReleaseImage(eigImage)
        cv.cvReleaseImage(tempImage)

    def __FindHarris(self, filename): #find the corners of images, and save all corner points in self.vKeyPoints
        self.img = highgui.cvLoadImage (filename)
        greyimg = cv.cvCreateImage(cv.cvSize(self.img.width, self.img.height), 8,1)
        w = cv.cvGetSize(self.img).width
        h = cv.cvGetSize(self.img).height
        
        image = cv.cvCreateImage(cv.cvGetSize(self.img), cv.IPL_DEPTH_32F, 1)
        cv.cvConvert(image, greyimg)
        self.cornerimg = cv.cvCreateImage(cv.cvGetSize(self.img), cv.IPL_DEPTH_32F, 1)
        cv.cvCornerHarris(image, self.cornerimg, 11,5,0.1)
        
    def SetBinary(self, t):
        self.drawimg = cv.cvCreateImage(cv.cvGetSize(self.img), 8, 3)
        cv.cvThreshold(self.img, self.drawimg, t, 255, cv.CV_THRESH_BINARY)
        
    def DrawKeyPoints(self):
        if (not self.drawimg):
            self.drawimg = cv.cvCloneImage(self.img) 

        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        ic = 0
        for c in self.points2:
            if (self.bDrawNumber):
                cv.cvPutText(self.drawimg, str(ic), cv.cvPoint(int(c.x), int(c.y)), myfont, cv.cvScalar(255, 255, 0,0))
                ic += 1
            cv.cvDrawCircle(self.drawimg, c, 4, cv.cvScalar(255,255,0,0))
        
    def SaveImage(self, filename):
        cv.cvNot(self.drawimg, self.drawimg)
        highgui.cvSaveImage(filename, self.drawimg)
 


    def GetCorner(self, fname, pointnum):
        self.nPointNum = pointnum
        #self.__FindHarris(fname)
        self.__FindCorner(fname)
    
    def GetCoordinate(self):
        X = []
        Y = []
        for k in self.vKeyPoints:
            X.append(k.x)
            Y.append(k.y)
        return X,Y    
        
              
    def __init__(self):
        self.output  = None
        self.bDraw  = 0
        self.bFilter = 1
        self.bDrawNumber = False
        self.drawimg = None
        
    
    
def main():
    ct = Corner()
    usage = "%s [options] <imgfile>" % (sys.argv[0])
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--display', action="store_true", dest = 'display', default = False, help = 'display the image')
    oparser.add_option('-n', '--number', dest = 'num', type='int',default = 200 , help = 'the number of feature points')
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')
    oparser.add_option('-s', '--save', dest = 'save', default = None, help = 'save the img file')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 2:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    ct.GetCorner(args[1], options.num)
    if (options.display):
        ct.DrawKeyPoints()
        highgui.cvNamedWindow ("Corner1", 1)
        highgui.cvShowImage ("Corner1", ct.drawimg)
        highgui.cvWaitKey (0)   
    if (options.save):
        highgui.cvSaveImage(options.save, ct.drawimg)    

if __name__ == '__main__':
    # load the image gived on the command line
    main()
