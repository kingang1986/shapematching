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
from array import * 
import math
import getopt
from contour import *
from TPS import *
import fileinput
from time import *
from aff import * 
from DPMatch import *
from aff import *



_green = cv.cvScalar (0, 255, 0, 0);
_blue = cv.cvScalar (255, 0, 0, 0);
_red = cv.cvScalar (0, 0, 255, 0);
_green1 = cv.cvScalar (128, 128, 0, 0);
_blue1 = cv.cvScalar (128, 0, 128, 0);
_red1 = cv.cvScalar (0, 128, 128, 0);
_green2 = cv.cvScalar (0, 255, 255, 0);
_blue2 = cv.cvScalar (255, 0, 255, 0);
_red2 = cv.cvScalar (0, 255, 255, 0);
_green3 = cv.cvScalar (128, 255, 255, 0);
_blue3 = cv.cvScalar (255, 128, 255, 0);
_red3 = cv.cvScalar (255, 255, 128, 0);

_black = cv.cvRealScalar (0)
_white = cv.cvRealScalar (255)

clrs = [_green, _blue, _red, _green1, _blue1, _red1, _green2, _blue2, _red2, _green3, _blue3, _red3 ] 

class SWMatch:

    def SubstituteCost(self, index1, index2):
        if (self.mFlag[index1][index2] == 1 ):
            return -10
        E1 = self.seq1[index1]
        E2 = self.seq2[index2]
        fCost1 = abs(sin(E1.angle/2 - math.pi/2) -  sin(E2.angle/2 - math.pi/2))  
        fCost2 = abs(E1.d2s - E2.d2s) / (E1.d2s + E2.d2s + 0.00001)
        fCost3 = abs(E1.d2c - E2.d2c) / (E1.d2c + E2.d2c + 0.00001)
        fCost31 = abs(E1.d2left - E2.d2left) / (E1.d2left + E2.d2left + 0.00001)
        fCost32 = abs(E1.d2right - E2.d2right) / (E1.d2right + E2.d2right + 0.00001)
        fCost4 = abs(E1.anglem/2 -  E2.anglem/2) / (abs(E1.anglem/2 - math.pi) +   abs(E1.anglem/2 - math.pi))
        fCost5 = abs(sin(E1.angle11 - math.pi) -  sin(E2.angle11 - math.pi))  
        fCost5 = abs(E1.angle11 - E2.angle11) / ( abs(E1.angle11 -math.pi) + abs(E2.angle11 - math.pi) + 0.001)  
        fCost6 = abs(E1.angle21/2 -  E2.angle21/2)  
        fCost7 = abs(E1.angle12/2 -  E2.angle12/2)  
        fCost8 = abs(E1.angle22/2 -  E2.angle22/2) 
        fCost9 = abs(E1.thickness -  E2.thickness) / ( E1.thickness + E2.thickness + 0.00001) 
        #print fCost1, fCost4, fCost2
        #return   0.2 - fCost1 - (fCost4 + (fCost2 * 5 ) + (fCost8 + fCost5 + fCost6 + fCost7))
        #print fCost5,
        #return min(0.3 - fCost2,  0.3 - fCost31, 0.3 - fCost32)
        w = max(E1.weight + 0.2, E2.weight + 0.2)  
        w = min(w, 1.0)
        w = max(w, 0.0)
        L2 = min(E1.d2s/self.avgd2s1, E2.d2s/self.avgd2s1)
        w2 = 0.5 / L2
        w2 = min(w2, 1.0)
        w3 = w + w2
        w = w / w3
        w2 = w / w3
        
        #return  (0.25 - fCost2) * (1 - w) +  (0.25 - fCost1) * w
        return self.costlim   -    fCost2  
 
    def GetDeleteCost(self, i, j1):
        return 0.3 
    def GetInsertCost(self, i1, j):
        return 0.3 
        
        
    def PrintDist(self):
        for L in self.mDist:
           for k in L:
               print "%3.3f "%(k),
           print ""
    def PrintPath(self):
        for L in self.mPath:
           for k in L:
               print "%3.3f "%(k),
           print ""              
    
    def InitImgFile(self, imgfile1, imgfile2):
        ct = Contour()
        ct.GetContour(imgfile1, 20)
        self.seq1 = ct.vKeyPoints[:]
        ct.DrawKeyPoints()
        self.img1 = ct.drawimg
        ct2 = Contour()
        ct2.GetContour(imgfile2, 20)
        if (self.bFlip):
          ct2.Flip()
        self.seq2 = ct2.vKeyPoints[:]
        ct2.DrawKeyPoints()
        self.img2 = ct2.drawimg
        self.GetXY()
        self.Flag1 = [0] * (self.m + 1)
        self.Flag2 = [0] * (self.n + 1)
        
        
        
    def InitFile(self, datafile1, datafile2):
        ct1 = Contour()
        ct1.Load(datafile1)
        self.seq1 = ct1.vKeyPoints[:]
        ct2 = Contour()
        ct2.Load(datafile2)
        if (self.bFlip):
          ct2.Flip()
        self.seq2 = ct2.vKeyPoints[:]
        self.GetXY()
        self.Flag1 = [0] * self.m
        self.Flag2 = [0] * self.n
        
          
    def InitData(self, X1, Y1):
        self.seq1 = []
        self.seq2 = []
        linecount = 0
        self.seq1 = X1[:]
        self.seq2 = Y1[:]
        self.m = len(self.seq1)
        self.n = len(self.seq2)
        self.GetXY()
        self.ExtendSeq1()

    def GetXY(self):
        self.m = len(self.seq1)
        self.n = len(self.seq2)

        self.X1 = []
        self.Y1 = []
        self.X2 = []
        self.Y2 = []
        self.avgd2s1 = 0
        self.avgd2s2 = 0
      
        
        for t in self.seq1:
            self.X1.append(t.x)
            self.Y1.append(t.y)
            self.avgd2s1 += t.d2s
            

        for t in self.seq2:
            self.X2.append(t.x)
            self.Y2.append(t.y)
            self.avgd2s2 += t.d2s
        self.avgd2s1 /= self.m
        self.avgd2s2 /= self.n

    def InitFlag(self):
        self.mFlag = []
        for i in range(self.m + 1):
            t = [0] * (self.n + 1)
            self.mFlag.append(t)
        
    def InitCostMatrix(self):
        self.mDist = []
        self.mPath = []
        self.mMoreCost = []
        for i in range(self.m + 1):
            d = array('f', [0.0] * (self.n + 1))
            self.mDist.append(d)
            k = array('i', [0] * (self.n + 1))
            self.mPath.append(k)
            s = array('i', [0] * (self.n + 1))
            self.mMoreCost.append(s)

    def CreateCostMatrix(self, start_r = 1, start_c = 1):
        m = self.m
        n = self.n
        self.InitCostMatrix()
        for i in range(start_r, (m + 1)):
            for j in range(start_c, (n + 1)):
                i1 = i - 1;
                j1 = j - 1;
                fSubstitute = self.mDist[i1][j1] + self.SubstituteCost(i1, j1) 
                fDelete = self.mDist[i1][j] - self.GetDeleteCost(i, j)
                fInsert = self.mDist[i][j1] - self.GetInsertCost(i, j)
                fNew = max([fSubstitute, fDelete, fInsert])
                #print  self.SubstituteCost(i1, j1),self.GetDeleteCost(i1, j1)/2,self.GetInsertCost(i1, j1)/2
                if (fNew > 0):
                    if (fSubstitute == fNew ): #/substitution is the best 
                        self.mPath[i][j] = 1
                        self.mMoreCost[i][j] = 0
                    elif (fDelete == fNew):
                        self.mPath[i][j] = 2
                        self.mMoreCost[i][j] = self.mMoreCost[i1][j] + 1
                    else:
                        self.mPath[i][j] = 3 
                        self.mMoreCost[i][j] = self.mMoreCost[i][j1] + 1
                fNew = max([fNew, 0])
                self.mDist[i][j] = fNew
   
        
    def GetSize(self):
        maxdist = -1
        #print len(self.X2), len(self.Y2)
        for i in range(len(self.X2)):
            for j in range(len(self.X2)):
                dist = (self.X2[i] - self.X2[j])  * (self.X2[i] - self.X2[j]) + (self.Y2[i] - self.Y2[j]) * (self.Y2[i] - self.Y2[j])
                if (dist > maxdist):
                    maxdist = dist
        self.size = math.sqrt(maxdist)          
        #print "size = %d, %f"%(len(self.X2), self.size) 
   
    def MatchNSeq(self, datafile1, datafile2):
        if (self.bImgfile):
            self.InitImgFile(datafile1, datafile2)
        else:
            self.InitFile(datafile1, datafile2)
        m = self.m;
        n = self.n;
        
        self.InitFlag()
        self.FindAnchor()
        self.CreateCostMatrix();
        fmax = 1
        fsum = 0
        self.nMatchPoint = 0
        #print self.m1, self.m, self.n
        while fmax > self.costlim:
            fmax = self.FindMatch()
            #print fmax
            if (fmax > self.costlim):
                fsum += fmax
            #self.PrintDist()
        #print "match points %d"%self.nMatchPoint
        #self.Show()
        return [fsum]

    def FindMatch(self):
        T = []
        for L in self.mDist:
            lmax = max(L)
            T.append(lmax)
        maxvalue = max(T)
        nRow = T.index(maxvalue)
        nCol = self.mDist[nRow].index(maxvalue)
        nMaxRow = nRow
        nMaxCol = nCol
        sub_cost = []
        X = []
        Y = []
        while (nRow > 0  or nCol > 0 ):
            t = self.mPath[nRow][nCol]
            if (t == 1) :
                f = self.mDist[nRow-1][nCol - 1]
                sub_cost.append(f)
                X.append(nRow - 1)
                Y.append(nCol - 1)
                nRow -= 1
                nCol -= 1
            elif (t == 2):
                nRow -= 1
            elif (t == 3):
                nCol -= 1
            elif (t == 0):
                break;
        if (self.bPrintPath):
            for i in range(len(X)):
                print "(%d -> %d) %f"%(X[i], Y[i], sub_cost[i])
        if (maxvalue > self.costlim * 1.5   and self.bDraw):
            mx1 = []
            my1 = []
            mx2 = []
            my2 = []
            for i in range(len(X)):
                mx1.append(self.X1[X[i]])
                my1.append(self.Y1[X[i]])
                mx2.append(self.X2[Y[i]])
                my2.append(self.Y2[Y[i]])
            self.DrawImg(mx1, my1, mx2, my2)        
        self.nMatchPoint += nMaxRow - nRow
        
        
        for i in range(nRow, self.m + 1):
            for j in range(0, nMaxCol+1):
                self.mFlag[i][j] = 1
        for i in range(0, nMaxRow + 1):
            for j in range(nCol, self.n + 1):
                self.mFlag[i][j] = 1
        self.CreateCostMatrix()
        
        return maxvalue
        
         
         
    def Chamfer(self):
        if (len(self.vMatchX) < 10):
            return [100000, 100000, 100000, 100000]
        PX1 =[]
        PY1 =[]
        PX2 =[]
        PY2 =[]
        for i in self.vMatchX:
           t = i % len(self.X1)
           PX1.append(self.X1[t])
           PY1.append(self.Y1[t])
        for i in self.vMatchY:
           PX2.append(self.X2[i])
           PY2.append(self.Y2[i])
        #print len(self.X1), len(self.X2)
        tp = TPS()
        tp.drawimg = self.bDraw
        bendingcost = tp.GetTransform(PX1, PY1, PX2, PY2)
        #DrawMatch(PX1, PY1, PX2, PY2)
        #print Afft
        chamfer = 0
        resultx = []
        resulty = []
        for i in range(len(self.X1)):
            ts = tp.Transform(self.X1[i], self.Y1[i])
            resultx.append(ts[0])
            resulty.append(ts[1])
            mindist = 1e10
            for j in range(len(self.X2)):
                dist =  math.sqrt((ts[0] - self.X2[j]) * (ts[0] - self.X2[j]) + (ts[1] - self.Y2[j]) * ( ts[1] - self.Y2[j]))
                if (dist < mindist):
                    mindist = dist
            chamfer += mindist
        #DrawMatch(self.X2, self.Y2, resultx, resulty)
        if (self.DrawImg):
          DrawIndexMatch(self.X1, self.Y1, self.X2, self.Y2, self.vMatchX, self.vMatchY)
        return bendingcost + [chamfer]
    def TPSMatch(self):
        mx1 = []
        my1 = []
        mx2 = []
        my2 = []
        for i in range(len(self.vMatchX)):
            mx1.append(self.X1[self.vMatchX[i]])
            my1.append(self.Y1[self.vMatchX[i]])
            mx2.append(self.X2[self.vMatchY[i]])
            my2.append(self.Y2[self.vMatchY[i]])
        tp = TPS()
        tp.drawimg = 0
        self.vTPS = tp.GetTransform(mx1, my1, mx2, my2)
        if self.bDraw:
            self.DrawImg(mx1, my1, mx2, my2)
        return self.vTPS
    def Show(self):
        highgui.cvNamedWindow ("image1", highgui.CV_WINDOW_AUTOSIZE)
        highgui.cvShowImage ("image1", self.img1)
        highgui.cvNamedWindow ("image2", highgui.CV_WINDOW_AUTOSIZE)
        highgui.cvShowImage ("image2", self.img2)
        highgui.cvWaitKey (0)
        #cv.cvNot(self.img1, self.img1)
        #cv.cvNot(self.img2, self.img2)
        highgui.cvSaveImage("match1.jpg", self.img1)
        highgui.cvSaveImage("match2.jpg", self.img2)

    def MatchDP(self):
        dp = DPMatch()
        #print self.Anchor1, self.Anchor2
        dp.seq1 = self.seq1[:self.Anchor1] + self.seq1[self.Anchor1:]
        dp.seq2 = self.seq2[:self.Anchor2] + self.seq2[self.Anchor2:]
        return dp.Match()
    def AffTrans(self):
        mx1 = []
        my1 = []
        mx2 = []
        my2 = []
        for i in range(len(self.vMatchX)):
            mx1.append(self.X1[self.vMatchX[i]])
            my1.append(self.Y1[self.vMatchX[i]])
            mx2.append(self.X2[self.vMatchY[i]])
            my2.append(self.Y2[self.vMatchY[i]])
        aff = AffTrans()
        re = aff.AffineTrans(mx1,my1,mx2,my2)
        dist = aff.GetDistance(mx1, my1, mx2, my2)
        if self.bDraw:
            self.DrawImg(mx1, my1, mx2, my2)

        return  re + [ dist] 
    def FindAnchor(self):
        m = self.m;
        n = self.n;
        self.InitCostMatrix()
        self.CreateCostMatrix();
        T = []
        for L in self.mDist:
            lmax = max(L)
            T.append(lmax)
        maxvalue = max(T)
        nRow = T.index(maxvalue)
        nCol = self.mDist[nRow].index(maxvalue)
        
        if (min(self.m - nRow , self.n -nCol) > 2):
            self.Anchor1 = nRow
            self.Anchor2 = nCol
            self.seq1 = self.seq1[nRow:] + self.seq1[:nRow]
            self.seq2 = self.seq2[nCol:] + self.seq2[:nCol]
            self.GetXY()
        else:
            while (nRow > 0  or nCol > 0 ):
                t = self.mPath[nRow][nCol]
                if (t == 1) :
                    nRow -= 1
                    nCol -= 1
                elif (t == 2):
                    nRow -= 1
                elif (t == 3):
                    nCol -= 1
                elif (t == 0):
                    break;
            self.Anchor1 = nRow
            self.Anchor2 = nCol
            self.seq1 = self.seq1[nRow:] + self.seq1[:nRow]
            self.seq2 = self.seq2[nCol:] + self.seq2[:nCol]
            self.GetXY()
           
        
    def MatchSeq(self, datafile1, datafile2):
        if (self.bImgfile):
            self.InitImgFile(datafile1, datafile2)
        else:
            self.bDraw = 0
            self.InitFile(datafile1, datafile2)
        m = self.m;
        n = self.n;
        
        self.InitFlag()
        self.InitCostMatrix()
        self.CreateCostMatrix();
        T = []
        for L in self.mDist:
            lmax = max(L)
            T.append(lmax)
        maxvalue = max(T)
        nRow = T.index(maxvalue)
        nCol = self.mDist[nRow].index(maxvalue)
        sub_cost = []
        self.vMatchX = []
        self.vMatchY  = []
        self.Anchor1 = nRow;
        self.Anchor2 = nCol;
        while (nRow > 0  or nCol > 0 ):
            t = self.mPath[nRow][nCol]
            if (t == 1) :
                f = self.mDist[nRow - 1][nCol - 1]
                if (self.mDist[nRow][nCol] > f):
                    sub_cost.append(self.mDist[nRow][nCol] - f)
                    self.vMatchX.append(nRow - 1)
                    self.vMatchY.append(nCol - 1)
                nRow -= 1
                nCol -= 1
            elif (t == 2):
                nRow -= 1
            elif (t == 3):
                nCol -= 1
            elif (t == 0):
                break;
        #for i in range(len(self.vMatchX)):
        #    print "(%d (%d %d) -> %d(%d %d)) %f"%(self.vMatchX[i], self.X1[self.vMatchX[i]],  self.Y1[self.vMatchX[i]], self.vMatchY[i],  self.X2[self.vMatchY[i]],  self.Y2[self.vMatchY[i]],sub_cost[i])
        self.Anchor1 += nRow
        self.Anchor2 += nCol
        self.Anchor1 /= 2
        self.Anchor2 /= 2
        result = [maxvalue]
        if (self.bDP):
            re = self.MatchDP()
            result = result + [re]
        if (self.bTPS):
            self.TPSMatch()
            result = result +  self.vTPS
        if (self.bAffTrans):
            result = result + self.AffTrans()
        if (maxvalue > 0.1 and self.bDraw):
            mx1 = []
            my1 = []
            mx2 = []
            my2 = []
            for i in range(len(self.vMatchX)):
                mx1.append(self.X1[self.vMatchX[i]])
                my1.append(self.Y1[self.vMatchX[i]])
                mx2.append(self.X2[self.vMatchY[i]])
                my2.append(self.Y2[self.vMatchY[i]])
            self.DrawImg(mx1, my1, mx2, my2)
            self.Show()        

        return result
        
    def DrawImg(self, X1, Y1, X2, Y2):
        for i in range(len(X1)):
            X1[i] = int(X1[i])
            Y1[i] = int(Y1[i])
        for i in range(len(X2)):
            X2[i] = int(X2[i])
            Y2[i] = int(Y2[i])
        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        for i in range(len(X1)):
            #cv.cvRectangle(img1, cv.cvPoint(X1[i] -1, Y1[i]-1), cv.cvPoint(X1[i]+1, Y1[i]+1), cv.cvScalar(255,0,255,0))
            cv.cvDrawCircle(self.img1, cv.cvPoint(X1[i], Y1[i]), 4, clrs[self.colorindex])
            cv.cvDrawCircle(self.img1, cv.cvPoint(X1[i], Y1[i]), 3, clrs[self.colorindex])
            cv.cvDrawCircle(self.img1, cv.cvPoint(X1[i], Y1[i]), 2, clrs[self.colorindex])
            cv.cvPutText(self.img1, str(i), cv.cvPoint(X1[i] + 8, Y1[i]), myfont, cv.cvScalar(255, 255,255, 0))
        for i in range(len(X2)):
            #cv.cvRectangle(img2, cv.cvPoint(X2[i]-1, Y2[i]-1), cv.cvPoint(X2[i]+1, Y2[i]+1), cv.cvScalar(255,0,255,0))
            cv.cvDrawCircle(self.img2, cv.cvPoint(X2[i], Y2[i]), 4, clrs[self.colorindex])
            cv.cvDrawCircle(self.img2, cv.cvPoint(X2[i], Y2[i]), 3, clrs[self.colorindex])
            cv.cvDrawCircle(self.img2, cv.cvPoint(X2[i], Y2[i]), 2, clrs[self.colorindex])
            cv.cvPutText(self.img2, str(i), cv.cvPoint(X2[i] + 8, Y2[i]), myfont, cv.cvScalar(255, 255, 255,0))
        self.colorindex += 1
        self.colorindex = self.colorindex % len(clrs)

    
    def __init__(self):
        self.m = 0
        self.n = 0
        self.vMatchX =[]
        self.vMatchY = []
        self.iSize = 1
        self.bDraw = 1
        self.bTPS = 1
        self.bFlip = 0  # flip image 1
        self.bImgfile = 0
        self.bPrintPath = 0
        self.normd2s = 0
        self.bDP = 1
        self.bAffTrans = 1
        self.colorindex = 0
        self.costlim = 0.25
        #self.drawimg = 0
        
def usage():
    print "usage: %s [-h] [-o] [-v] [-d] [-n] [-i] [-f] datafile1 datafile2"%sys.argv[0]
    
def main():
   
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:vdtif", ["help", "output=", "draw", "tps", "image", "flip"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    output = None
    matcher = SWMatch()
    verbose = False
    bDraw = 1
    npoint = 20
    bTPS = 1
    for o, a in opts:
        if o == "-v":
           verbose = 1
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--output"):
            output = a
        if o in ("-d", "--draw"):
            bDraw = 0
        if o in ("-t", "--tps"):
            bTPS = 0
        if o in ("-i", "--image"):
            matcher.bImgfile = 1
        if o in ("-f", "--flip"):
            matcher.bFlip = 1
    if (len(args)) != 2:
        usage()
        sys.exit(2)

    
    matcher.bDraw = bDraw
    matcher.bTPS = 0
    cost = matcher.MatchSeq(args[0], args[1])
    
    #matcher.PrintDist()
    print cost  
    #print matcher.vTPS
    #highgui.cvSaveImage("match_tmp.bmp", matcher.drawimg)
  
if __name__ == '__main__':
    # load the image gived on the command line
    main()
    
     
    
