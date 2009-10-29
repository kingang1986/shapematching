#!/usr/bin/env python
"""MSS_SWMatch - a Python class
"""
__author__ = "Longbin Chen(longbinc@yahoo.com)"
__version__ = "$Revision: 2.0 $"
__history__="""
    Created by Longbin Chen
    Modified 08/30/2009, by longbinc@yahoo.com
         1. change file format
    """

import sys, os, string, optparse, fileinput, string, math, numpy
from MSS import *
from exMSS import *
from exDAS import * 
from exAngle import *
from exSC import *


OP_MATCH = 1
OP_SUBST = 2
OP_DELET = 3 
OP_INSRT = 4 
NO_OP =  -1
NO_OP_BOTH = -1
NO_OP_A = -2
NO_OP_B = -3 


def PlotMatch(bX, bY, matchscore, file1, file2):
    x1,y1 = loadxy(file1)
    x2,y2 = loadxy(file2)
    for i in range(len(x2)):
        x2[i] += 800
    printheader(len(x1))
    m1 = []
    for i in range(len(x1)):
        m1.append(1)
    m2 = []
    for i in range(len(x2)):
        m2.append(1)
    for i in range(len(bX)):
        if (bY[i] == -1 and bX[i] != -1):
            m1[bX[i]] = 0
    for i in range(len(bY)):
        if (bX[i] == -1 and bY[i] != -1):
            m2[bY[i]] = 0 
    wei1 = []
    wei2 = []
    for i in range(len(x1)):
       wei1.append(8)
       wei2.append(8)
    sumsc = 0
    cnt = 0
    for t in matchscore:
        if t > 0:
            sumsc += t
            cnt +=1 
    sumsc = sumsc/cnt
    normscore = [ss / sumsc for ss in matchscore]
    for i in range(len(normscore)):
       if (normscore[i] < 0.5):
           sc = 8
       else:
           sc = 8  + (normscore[i] - 0.5) * 20
       wei1[bX[i]] = sc
       wei2[bY[i]] = sc
    print >> sys.stderr, len(x1), len(y1), len(matchscore), len(normscore), len(wei1), len(wei2)
    printWeightPts(x1, y1, 'n', m1, wei1)
    printWeightPts(x2, y2, 'm', m2, wei2)
    edx = []
    edy = []
    for i in range(len(bX)):
        if (bX[i] != -1 and bY[i] != -1):
            edx.append(bX[i])
            edy.append(bY[i])   
    printEdge(edx, edy)
    printtail()


class SmithWaterman:
    
    def __init__(self):
        self.weights = None
        self.matchmap = None
        self.gapmap = None
        self.feature_length = 0

    def SetModel(self, weights, matchmap, gapmap):
        self.weights = weights[:]
        self.matchmap = matchmap[:]
        self.gapmap = gapmap[:]
        
    def GetSubScore(self, index1, index2):
    	f = 0.0
        fkeys = self.seq1.points[0].getvalues().keys()
        d = 0 
        if (self.matchmap == None):
            for k in fkeys: 
    	        d  =  abs(self.seq1.points[index1].getvalue(k) -  self.seq2.points[index2].getvalue(k))  
                m  =  abs(self.seq1.points[index1].getvalue(k)) +  abs(self.seq2.points[index2].getvalue(k)) 
                if ( m > 0):
                    f += d  / m 
        f = 20 - f
        return f

    def GetGapScore(self, data1):
        return  1.0
 
    def PrintPath(self, m):
        for L in m:
           for k in L:
               print "%3.3f "%(k),
           print ""              
    def Align(self, seq1, seq2):
    	self.seq1 = seq1
    	self.seq2 = seq2
    	self.feature_length = len([f for f in seq1.points[0].getvalues() if f.startswith("f")])
    	MatchX = []
    	MatchY = []
        matchscores = []
        m = len(seq1.points) 
        n = len(seq2.points)
        cost = numpy.zeros((m + 1, n + 1))
        path = numpy.zeros((m + 1, n + 1))
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                i1 = i -1
                j1 = j -1
                fSubstitute = cost[i1][j1] + self.GetSubScore(i1, j1) 
                fDelete = cost[i1][j] + self.GetGapScore(i1)
                fInsert = cost[i][j1] + self.GetGapScore(j1)
                fNew = max([fSubstitute, fDelete, fInsert])
                #print  self.SubstituteCost(i1, j1),self.GetDeleteCost(i1, j1)/2,self.GetInsertCost(i1, j1)/2
                if (fNew > 0):
                    if (fSubstitute == fNew ): #/substitution is the best 
                        path[i][j] = 1
                    elif (fDelete == fNew):
                        path[i][j] = 2
                    else:
                        path[i][j] = 3 
                fNew = max([fNew, 0])
                cost[i][j] = fNew  
        aligncost = cost.max()
        npos = numpy.argmax(cost)
        nRow = npos / (n + 1)
        nCol = npos % (n + 1)
        maxrow = nRow
        maxcol = nCol
        alignMatch = []
        subsvalue = []
        while (nRow > 0 and nCol > 0 ):
            t = path[nRow][nCol]
            if (t == 1):
                MatchX.append(nRow - 1)
                MatchY.append(nCol - 1)
                matchscores.append(cost[nRow][nCol] - cost[nRow-1][nCol-1])
                nRow -= 1
                nCol -= 1
                alignMatch.append(OP_SUBST)
            elif (t == 2):
                MatchX.append(nRow - 1)
                MatchY.append(- 1)
                matchscores.append(-1)
                nRow -= 1
                alignMatch.append(OP_DELET)
            elif (t == 3):
                MatchX.append(- 1)
                MatchY.append(nCol - 1)
                matchscores.append(-1)
                nCol -= 1
                alignMatch.append(OP_INSRT)
            elif (t == 0):
                break
        return aligncost, alignMatch, MatchX, MatchY, matchscores


class MSS_SW:
    
    def align(self, mss1, mss2):
        sumscore = []
        matcher = SmithWaterman()
        idx = -1
        bXFinal = []
        bYFinal = []
        bScoreFinal = []
        mC = []
        for i1, c1 in enumerate(mss1.seqs):
            idx += 1
            cscore = -100000000
            bestcurve = None
            for i2, c2 in enumerate(mss2.seqs):
                cost, align, X, Y, mscores = matcher.Align(c1, c2)
                print X, Y
                if (cost > cscore):
                    #print >>sys.stderr, cost, align, len(X), len(Y), len(mscores)
                    cscore = cost
                    #print >>sys.stderr, "---------------"
                    bScore = mscores[:]
                    bX = []
                    bY = []
                    for x in X:
                       if x != -1: 
                          bX.append((i1, x))
                       else:
                          bX.append((i1, -1))
                    for y in Y:
                       if y != -1:
                          bY.append((i2, y))
                       else:
                          bY.append((i2, -1))
                    bestcurve = c2
            sumscore.append(cscore)   
            bXFinal += bX
            bYFinal += bY
            bScoreFinal += bScore
        return sumscore, bXFinal, bYFinal,bScoreFinal
 

def main():

    usage = "%prog [options] <mss_file1> <mss_file2>"
    version = "%prog 0.2\nLongbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')
    oparser.add_option('-i', '--image', dest = 'image', action="store_true", default = False, help = 'input image file')
    oparser.add_option('-m', '--drawnumber', action="store_true", dest = 'drawnumber', default = False, help = 'display the point numbers')
    oparser.add_option('-t', '--threshold', dest = 'threshold', type='int',default = 100 , help = 'the threshold for image binarification')
    oparser.add_option('-n', '--number', dest = 'num', type='int',default = 100 , help = 'the number of feature points')


    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    if (options.image):
        das = exDAS()
        agl = exAngle()
        sc = exSC()
        ct1 = ExtractMSS()
        ct1.GetContour(args[1], options)
        ct1.DrawKeyPoints()
        mss1 = ct1.mss
        for seq in mss1.seqs:
             das.ExtractFeature(seq.points, None)
             agl.ExtractFeature(seq.points, None)
        sc.ExtractFeature(mss1)
        
        mss1.save("_1.tmp")

        ct2 = ExtractMSS()
        ct2.GetContour(args[2], options)
        ct2.DrawKeyPoints()
        mss2 = ct2.mss
        for seq in mss2.seqs:
             das.ExtractFeature(seq.points, None)
             agl.ExtractFeature(seq.points, None)
        sc.ExtractFeature(mss2)
        mss1.save("_2.tmp")
    else:
        mss1 = MSS()
        mss1.load(args[1])
        mss2 = MSS()
        mss2.load(args[2])
    sw = MSS_SW()
    sumsc, bX, bY, sfinal =  sw.align(mss1, mss2)
    print sumsc[0]
    if (options.image):
        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        ptcount = 0
        for i in range(len(bX)):
            xs, xi = bX[i]
            ys, yi = bY[i]
            if (xi != -1 and yi != -1):
                    ptcount  += 1
 #                   cv.cvDrawCircle(ct1.drawimg, cv.cvPoint(int(mss1.seq[xs].points[xi].x), int(mss1.seq[xs].points[xi].y)),2, clrs[idx])
                    cv.cvPutText(ct1.drawimg, str(ptcount), cv.cvPoint(int(mss1.seqs[xs].points[xi].x), int(mss1.seqs[xs].points[xi].y)), myfont,  cv.cvScalar(255,255,255,0))
#                    cv.cvDrawCircle(ct2.drawimg, cv.cvPoint(int(bestcurve[yi].x), int(bestcurve[yi].y)),2, clrs[idx])
                    cv.cvPutText(ct2.drawimg, str(ptcount), cv.cvPoint(int(mss2.seqs[ys].points[yi].x), int(mss2.seqs[ys].points[yi].y)), myfont,  cv.cvScalar(255,255,255,0))
        highgui.cvNamedWindow ("contour1", 1)
        highgui.cvNamedWindow ("contour2", 1)
        highgui.cvShowImage ("contour1", ct1.drawimg)
        highgui.cvShowImage ("contour2", ct2.drawimg)
        highgui.cvWaitKey (0)
    if options.output:
        print bX, bY
        fout = open(options.output, 'w')
        for it, b in enumerate(bX):
            fout.write("%d\t%d\t%d\t%d\n" %(b[0], b[1], bY[it][0],bY[it][1]))
        fout.close()
        
        
if __name__ == '__main__':
    # load the image gived on the command line
    main()
