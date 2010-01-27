#! /usr/bin/env python
'''
extract angle features from the curves, 
7 features are extracted:
fag11, fag22, fag33, fag44 : angle with different granulatrity

 
'''

import sys, math, string, fileinput, optparse
from MSS import *

class exAngle:
    def __GetB(self, a, b, c):
        xa = self.keypoints[a].x 
        xb = self.keypoints[b].x 
        xc = self.keypoints[c].x 
        ya = self.keypoints[a].y 
        yb = self.keypoints[b].y 
        yc = self.keypoints[c].y 
        x0 = (xa + xc) /2.0
        y0 = (ya + yc) /2.0
        xa -= x0
        xb -= x0
        xc -= x0
        ya -= y0
        yb -= y0
        yc -= y0
        L = math.sqrt( xc * xc + yc * yc)
        x = xb * xc / L + yb * yc / L 
        y = - xb * yc / L  + yb * xc / L  
        return x, y

    def GetBlur(self, distn ):
        nPoint = len(self.keypoints)
        for i in range(nPoint):
            self.keypoints[i].fbx = 0
            self.keypoints[i].fby = 0
            self.keypoints[i].gby = 0
            self.keypoints[i].gbx = 0
        for i in range(distn, nPoint - distn):
            dx, dy = self.__GetB(i - distn, i, i + distn) 
            self.keypoints[i].fbx = dx
            self.keypoints[i].fby = dy
            self.keypoints[i].gbx = -dx
            self.keypoints[i].gby = -dy
        
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
            self.keypoints[i].gag11 = 0
            self.keypoints[i].gag22 = 0
            self.keypoints[i].gag33 = 0
            self.keypoints[i].gag44 = 0
            self.keypoints[i].gagout = 0
            self.keypoints[i].gagin = 0
            self.keypoints[i].gcurve = 0
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
            self.keypoints[i].gag11 = - self.keypoints[i].fag11
            self.keypoints[i].gag22 = - self.keypoints[i].fag22
            self.keypoints[i].gag33 = - self.keypoints[i].fag33
            self.keypoints[i].gag44 = - self.keypoints[i].fag44
            self.keypoints[i].gagout = max(0,  -self.keypoints[i].gag11)
            self.keypoints[i].gagin = max(0,  self.keypoints[i].gag11)
            
        for i in range(nPoint):
            self.keypoints[i].fcurve = (self.keypoints[(i + 2) % nPoint].fag22 + self.keypoints[ (i - 2 + nPoint) % nPoint].fag22 - 2 * self.keypoints[i].fag22) / 2.0
            self.keypoints[i].gcurve = -self.keypoints[i].fcurve

 
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
        self.GetBlur(2)
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

