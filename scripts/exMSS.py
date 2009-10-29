#!/usr/bin/env python

import sys, math, string, optparse, fileinput

from MSS import *
from opencv import cv
from opencv import highgui


class ExtractMSS:

    def __init__(self):
        self.output  = None
        self.bDraw  = 0
        self.bFilter = 1
        self.bDrawNumber = 1
        self.threshold = 100 
        self.mss = MSS()

    def __findContour(self, filename): #find the contour of images, and save all points in self.vKeyPoints
        self.img = highgui.cvLoadImage (filename)
        self.grayimg = cv.cvCreateImage(cv.cvSize(self.img.width, self.img.height), 8,1)
        self.drawimg = cv.cvCreateImage(cv.cvSize(self.img.width, self.img.height), 8,3)
        cv.cvCvtColor (self.img, self.grayimg, cv.CV_BGR2GRAY)
        cv.cvSmooth(self.grayimg, self.grayimg, cv.CV_BLUR, 9)
        cv.cvSmooth(self.grayimg, self.grayimg, cv.CV_BLUR, 9)
        cv.cvSmooth(self.grayimg, self.grayimg, cv.CV_BLUR, 9)
        cv.cvThreshold( self.grayimg, self.grayimg, self.threshold, self.threshold +100, cv.CV_THRESH_BINARY )
        cv.cvZero(self.drawimg)
        storage = cv.cvCreateMemStorage(0)
        nb_contours, cont = cv.cvFindContours (self.grayimg,
            storage,
            cv.sizeof_CvContour,
            cv.CV_RETR_LIST,
            cv.CV_CHAIN_APPROX_NONE,
            cv.cvPoint (0,0))
            
        cv.cvDrawContours (self.drawimg, cont, cv.cvScalar(255,255,255,0), cv.cvScalar(255,255,255,0), 1, 1, cv.CV_AA, cv.cvPoint (0, 0))
        self.allcurve = []
        idx = 0
        for c in cont.hrange():
            PointArray = cv.cvCreateMat(1, c.total  , cv.CV_32SC2)
            PointArray2D32f= cv.cvCreateMat( 1, c.total  , cv.CV_32FC2)
            cv.cvCvtSeqToArray(c, PointArray, cv.cvSlice(0, cv.CV_WHOLE_SEQ_END_INDEX))
            fpoints = []
            for i in range(c.total):
                kp = myPoint()
                kp.x = cv.cvGet2D(PointArray,0, i)[0]
                kp.y = cv.cvGet2D(PointArray,0, i)[1]
                kp.index = idx
                idx += 1
                fpoints.append(kp)
            self.allcurve.append(fpoints)
        self.curvelength = idx

    def __GetCurve(self, iOrder = 2):
        for c in self.allcurve:    
            nPoint = len(c)
            for i in range(nPoint):
                #print (i + iOrder) % nPoint, (i - iOrder) % nPoint, i, nPoint, len(self.vKeyPoints)
                f = (c[(i + iOrder) % nPoint].angle + c[ (i - iOrder + nPoint) % nPoint].angle - 2 * c[i].angle) / 2.0
                self.vKeyPoints[i].curve = f
            
    def __distance(self, x, y, nlength):
        sm = min(x, y) 
        larg = max(x, y) 
        L1 = larg - sm
        L2 = sm + nlength - larg
        return min(L1, L2) 

    def __evenSample(self, npoint):
        ndist = self.curvelength / (0.001 + npoint)
        self.mss = MSS()
        for c in self.allcurve:
            seq = Sequence()
            for p in c:
                i = p.index
                #print int((i + 1) / ndist), int(i/ ndist)
                if (int((i + 1) / ndist) - int(i / ndist)) == 1:
                    seq.points.append(p)
            self.mss.seqs.append(seq)

    def GetContour(self, fname, options):
          self.bDrawNumber = options.drawnumber
          self.threshold = options.threshold
          self.__findContour(fname)
          self.__evenSample(options.num)

    def SaveImage(self, filename):
        cv.cvNot(self.drawimg, self.drawimg)
        highgui.cvSaveImage(filename, self.drawimg)
 
        
    def DrawKeyPoints(self):
        ic = 0
        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        for ic, c in enumerate(self.mss.seqs):
          for k in c.points:
            if (self.bDrawNumber):
                cv.cvPutText(self.drawimg, str(ic), cv.cvPoint(int(k.x), int(k.y)), myfont, cv.cvScalar(255, 255, 0,0))
            cv.cvDrawCircle(self.drawimg, cv.cvPoint(int(k.x), int(k.y)), 4, cv.cvScalar(255,0,255,0))
            #cv.cvDrawCircle(self.drawimg, cv.cvPoint(int(k.x), int(k.y)), 6, cv.cvScalar(255,255,255,0))
            #cv.cvDrawCircle(self.drawimg, cv.cvPoint(int(k.x), int(k.y)), 5, cv.cvScalar(255,255,255,0))
            #cv.cvDrawCircle(self.drawimg, cv.cvPoint(int(k.x), int(k.y)), 2, cv.cvScalar(255,255,255,0))
    
    
def main():
  
    usage = "%prog [options] <imgfile>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--display', action="store_true", dest = 'display', default = False, help = 'display the image')
    oparser.add_option('-m', '--drawnumber', action="store_true", dest = 'drawnumber', default = False, help = 'display the point numbers')
    oparser.add_option('-n', '--number', dest = 'num', type='int',default = 200 , help = 'the number of feature points')
    oparser.add_option('-t', '--threshold', dest = 'threshold', type='int',default = 100 , help = 'the threshold for image binarification')
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')
    oparser.add_option('-s', '--save', dest = 'save', default = None, help = 'save the img file')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 2:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
        
    ct = ExtractMSS()
    ct.GetContour(args[1], options)

    if (options.display):
        ct.DrawKeyPoints()
        highgui.cvNamedWindow ("contour", 1)
        highgui.cvShowImage ("contour", ct.drawimg)
        highgui.cvWaitKey (0)       

    if (options.output):
        ct.mss.save(options.output)

    if (options.save):
        highgui.cvSaveImage(options.save, ct.drawimg)    

if __name__ == '__main__':
    main()
