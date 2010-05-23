#! /usr/bin/env python
 

# import the necessary things for OpenCV

import sys, math, string, fileinput, optparse
from opencv import cv
from opencv import highgui
from MSS import *
 


class exDAS:
    def __GetShapeCenter(self):
        self.fAvgX = 0.0
        self.fAvgY = 0.0
        npoint = len(self.keypoints)
        if npoint > 0:
	    for k in self.keypoints:
	        self.fAvgX += k.x
	        self.fAvgY += k.y
	    self.fAvgX /= npoint
	    self.fAvgY /= npoint
                
    def __GetAngleBetween(self, p0, p1, p2):
        dx1 = p0.x - p1.x
        dy1 = p0.y - p1.y
        dx2 = p2.x - p1.x
        dy2 = p2.y - p1.y
        a1 = math.atan2(dx1, dy1)
        a2 = math.atan2(dx2, dy2)
        da = a1 - a2
        if (da < 0):
           da += 2 * math.pi
        if (da > 2 * math.pi):
            da -= 2 * math.pi
        return da
    def __GetCentAngle(self, p0, p1, p2):
        dx1 = p0.x - p1.x
        dy1 = p0.y - p1.y
        dx2 = p2.x - p1.x
        dy2 = p2.y - p1.y
        a1 = math.atan2(dx1, dy1)
        a2 = math.atan2(dx2, dy2)
        da = a1 - a2
        if (da < 0):
           da += 2 * math.pi
        if (da > 2 * math.pi):
            da -= 2 * math.pi
        ca = a2 + da /2 
        if (ca < 0):
            ca += 2 * math.pi 

        if (ca > 2 * math.pi):
            ca -= 2 * math.pi
     
        return ca

    def __GetLeftAngle(self, p0, p1, p2):
        dx = p0.x - p1.x
        dy = p0.y - p1.y
        a1 = math.atan2(-dy, dx)
        return a1

    def __GetRightAngle(self, p0, p1, p2):
        dx = p1.x - p2.x
        dy = p1.y - p2.y
        a1 = math.atan2(-dy, dx)
        return a1

    def __GetDist(self, i, j):
    	dx = self.keypoints[i].x - self.keypoints[j].x
    	dy = self.keypoints[i].y - self.keypoints[j].y
    	return math.sqrt(dx * dx + dy * dy) / self.fScale

    def __GetSupport(self):
        nPoint = len(self.keypoints)
        qrt = nPoint /4
        for i in range(nPoint):
            i0 =  (i - qrt + nPoint) % nPoint
            i1 =  (i + qrt) % nPoint 
            i2 =  (i + qrt * 2) % nPoint 
            self.keypoints[i].fspt1 = self.__GetDist(i0, i)
            self.keypoints[i].fspt2 = self.__GetDist(i1, i)
            self.keypoints[i].fspt3 = self.__GetDist(i2, i)
 
    def __GetDDAS(self):
        nPoint = len(self.keypoints)
        for i in range(nPoint):
            i0 =  (i - 1 + nPoint) % nPoint
            i2 =  (i + 1) % nPoint 
            self.keypoints[i].fddleft1 = self.keypoints[i0].fd2left - self.keypoints[i].fd2left
            self.keypoints[i].fddleft2 = self.keypoints[i2].fd2left - self.keypoints[i].fd2left
            self.keypoints[i].fddright1 = self.keypoints[i0].fd2right - self.keypoints[i].fd2right
            self.keypoints[i].fddright2 = self.keypoints[i2].fd2right - self.keypoints[i].fd2right
            
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
        da = a1 - a2

        if (da < 0):
                da += math.pi * 2
        if (da >  2* math.pi ):
                da -= math.pi * 2
        return da
    def __IsClockwise(self):
       nPoint = len(self.keypoints)
       fsumangle = 0.0
       for i in range(nPoint):
            p0 = (i - 1 + nPoint) % nPoint
            p1 = i
            p2 = (i + 1) % nPoint
            fCentAngle = self.__Get3PAngle(p0, p1, p2)
            fsumangle += fCentAngle
       if abs(fsumangle - (nPoint -2 ) * 3.1415926) < 0.1:
           return 1
       else:
           return 0
      
    def __GetDAS(self):
        nPoint = len(self.keypoints)
        for i in range(nPoint):
            p0 = self.keypoints[(i - 1 + nPoint) % nPoint]
            p1 = self.keypoints[i]
            p2 = self.keypoints[(i + 1) % nPoint]
            fCentAngle = self.__GetCentAngle(p0, p1, p2)
            dx1 = math.sin(fCentAngle)
            dy1 = math.cos(fCentAngle)
            fBisectDist = self.__GetCrossDist(p1, dx1, dy1, i)
            self.keypoints[i].fdas = fBisectDist / self.fScale
            angle = self.__GetAngleBetween(p0,p1,p2)
            if (angle > math.pi / 2.0):
                fLeftAngle = self.__GetLeftAngle(p0, p1, p2)
                dx1 = math.sin(fLeftAngle)
                dy1 = math.cos(fLeftAngle)
                fLeftDist = self.__GetCrossDist(p1, dx1, dy1, i)
                self.keypoints[i].fd2left = fLeftDist / self.fScale
                fRightAngle = self.__GetRightAngle(p0, p1, p2)
                dx1 = math.sin(fRightAngle)
                dy1 = math.cos(fRightAngle)
                fRightDist = self.__GetCrossDist(p1, dx1, dy1, i)
                self.keypoints[i].fd2right = fRightDist /self.fScale
            else:
                self.keypoints[i].fd2left = 0
                self.keypoints[i].fd2right = 0
            dx = (self.fAvgX - p1.x)
            dy = (self.fAvgY - p1.y)
            self.keypoints[i].fd2c = math.sqrt(dx * dx + dy * dy) / self.fScale  

    def __GetCrossDist(self, p1, dx, dy, iPointIndex):
        bFound = 0
        fDist = 0
        bestPoint = cv.cvPoint(0, 0)
        bestLength = 1e10
        bigLength = -1
        nPoints = len(self.keypoints)
        for k in range(nPoints):
            if (k == iPointIndex or k == iPointIndex + 1):
                continue
            q1 = self.keypoints[(k - 1 + nPoints) % nPoints]
            q2 = self.keypoints[k]
            du = q2.x - q1.x
            dv = q2.y - q1.y
            dd = (dy * du - dx * dv)
            if (dd == 0):
                continue
            t =  (dy * (p1.x - q1.x) - dx * (p1.y - q1.y)) / dd
            if (t >= -0.0001 and t <= 1.0001): # found it
                ptt =  cv.cvPoint(int(q1.x + t * du), int(q1.y + t * dv))
                l = math.sqrt((ptt.x - p1.x ) * (ptt.x - p1.x ) + (ptt.y - p1.y ) * (ptt.y - p1.y))
                l2 = ((dv * q1.x - du * q1.y) - (dv * p1.x - du * p1.y)) / ( dv * dx - du * dy)
                bFound = 1
                if (l <= bestLength and l2 > 0):
                    bestPoint = ptt
                    bestLength = l
        fDist = bestLength
        if (not bFound):
            fDist = 0
        if (self.img):
            cv.cvLine(self.img, cv.cvPoint(int(p1.x), int(p1.y)), bestPoint, cv.cvScalar(255, 255, 255, 0))
        return fDist

    def __GetAverageDist(self):
        fdist = 0
        for j in self.keypoints:
            d = math.sqrt( (j.x - self.fAvgX) * (j.x - self.fAvgX) + (j.y - self.fAvgY) * (j.y - self.fAvgY))
            fdist += d
        fdist /= len(self.keypoints)
        if (fdist == 0):
            fdist = 1
        self.fScale = fdist
       
    
    def ExtractFeature(self, keypoints, img):
        self.keypoints = keypoints
        self.img = img
        self.__GetShapeCenter()
        self.__GetAverageDist()
        if (self.__IsClockwise()):
            self.__GetDAS()
        else:
            self.keypoints.reverse()
            self.__GetDAS()
            self.keypoints.reverse()
            for i in range(len(self.keypoints)):
                self.keypoints[i].fdas = - self.keypoints[i].fdas
                self.keypoints[i].fd2left = - self.keypoints[i].fd2left
                self.keypoints[i].fd2right = - self.keypoints[i].fd2right
        self.__GetDDAS()
        self.__GetSupport()
        return keypoints
         
def main():

    usage = "%prog [options] <msshape> <outputmss>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    #oparser.add_option('-c', '--shapecontext', action="store_false",dest = 'shapecontext', default = True,  help = 'extract shape context  features(default true)')
    #oparser.add_option('-w', '--das', action="store_false",dest = 'das', default = True,  help = 'extract DAS features(default true)')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)

    das = exDAS()
    mss = MSS()
    mss.load(args[1])
    for seq in mss.seqs:
        das.ExtractFeature(seq.points, None)
    mss.save(args[2])

if __name__ == '__main__':
    # load the image gived on the command line
    main()
 
