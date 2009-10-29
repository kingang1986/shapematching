#! /usr/bin/env python
 

# import the necessary things for OpenCV

import sys, math, string, fileinput, optparse
from opencv import cv
from opencv import highgui
from MSS import *
  

class exAngle:
    def __GetNAngle(self, distn):
        nPoint   = len(self.keypoints)
        for i in range(nPoint):
            self.keypoints[i].fag11 = 0
            self.keypoints[i].fag22 = 0
            self.keypoints[i].fag33 = 0
            self.keypoints[i].fag44 = 0
            self.keypoints[i].fagout = 0
            self.keypoints[i].fagin = 0
            self.keypoints[i].fcurve = 0
        for i in range(nPoint):
            iPrev = (i - distn + nPoint) % nPoint
            iNext = (i + distn + nPoint) % nPoint
            self.keypoints[i].fag11 = self.__Get3PAngle(iPrev, i, iNext)
            iPrev = (i - 2 * distn + nPoint) % nPoint
            iNext = (i + 2 * distn + nPoint) % nPoint
            self.keypoints[i].fag22 = self.__Get3PAngle(iPrev, i, iNext)
            iPrev = (i - 3 * distn + nPoint) % nPoint
            iNext = (i + 3 * distn + nPoint) % nPoint
            self.keypoints[i].fag33 = self.__Get3PAngle(iPrev, i, iNext)
            iPrev = (i - 4 * distn + nPoint) % nPoint
            iNext = (i + 4 * distn + nPoint) % nPoint
            self.keypoints[i].fag44 = self.__Get3PAngle(iPrev, i, iNext)
            self.keypoints[i].fagout = max(0, - self.keypoints[i].fag11)
            self.keypoints[i].fagin = max(0, self.keypoints[i].fag11)
        for i in range(nPoint):
            self.keypoints[i].fcurve = (self.keypoints[(i + 2) % nPoint].fag22 + self.keypoints[ (i - 2 + nPoint) % nPoint].fag22 - 2 * self.keypoints[i].fag22) / 2.0

 
    def __GetAngle(self):
        nPoint   = len(self.keypoints)
        for i in range (nPoint):
            iPrev = (i - 1 + nPoint) % nPoint 
            iNext = (i + 1 + nPoint) % nPoint 
            self.keypoints[i].fag11 = self.__Get3PAngle(iPrev, i, iNext)
            iPrev = (i - 2 + nPoint) % nPoint 
            iNext = (i + 2 + nPoint) % nPoint 
            self.keypoints[i].fag22 = self.__Get3PAngle(iPrev, i, iNext)
            iPrev = (i - 3 + nPoint) % nPoint 
            iNext = (i + 3 + nPoint) % nPoint 
            self.keypoints[i].fag33 = self.__Get3PAngle(iPrev, i, iNext)
            iPrev = (i - 4 + nPoint) % nPoint 
            iNext = (i + 4 + nPoint) % nPoint 
            self.keypoints[i].fag44 = self.__Get3PAngle(iPrev, i, iNext)
            self.keypoints[i].fagout = max(0, - self.keypoints[i].fag11)
            self.keypoints[i].fagin = max(0, self.keypoints[i].fag11)
        for i in range(nPoint):
            self.keypoints[i].fcurve = (self.keypoints[(i + 2) % nPoint].fag22 + self.keypoints[ (i - 2 + nPoint) % nPoint].fag22 - 2 * self.keypoints[i].fag22) / 2.0
            
    def __Get3PAngle(self, p0, p1, p2):
        x0 = self.keypoints[p0].x
        y0 = self.keypoints[p0].y
        x1 = self.keypoints[p1].x
        y1 = self.keypoints[p1].y
        x2 = self.keypoints[p2].x
        y2 = self.keypoints[p2].y
        dx1 = x0 - x1 
        dy1 = y0 - y1 
        dx2 = x2 - x1 
        dy2 = y2 - y1  
        a1 = math.atan2(dx1, dy1)
        a2 = math.atan2(dx2, dy2)
        da = a1 - a2 - math.pi
        
        if (da < -math.pi):
                da += math.pi * 2
        if (da >  math.pi ):
                da -= math.pi * 2
        return da

    def ExtractFeature(self, keypoints, img):
        self.keypoints = keypoints
        self.img = img
        self.__GetNAngle(3)
        return keypoints
        
        
    
def main():

    usage = "%prog [options] <msshape> <outputmss>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)

    agl = exAngle()
    mss = MSS()
    mss.load(args[1])
    for seq in mss.seqs:
        agl.ExtractFeature(seq.points, None)
    mss.save(args[2])

if __name__ == '__main__':
    # load the image gived on the command line
    main()

