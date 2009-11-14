#!/usr/bin/env python
#####################################################
# Longbin Chen
# ------------
# Created by longbinc@yahoo.com
#####################################################

import sys, math, string, optparse, fileinput, cStringIO, random

from opencv import cv
from opencv import highgui

from MSS import *
import NW

_red = cv.cvScalar (0, 0, 255, 0)
_green = cv.cvScalar (0, 255, 0, 0)
_white = cv.cvScalar (255,255,255,0)
_black = cv.cvScalar (0,0,0,0)
color = [_red, _green, _white, cv.cvScalar(128,0, 128), cv.cvScalar(128,128, 128), cv.cvScalar(128,0, 0), cv.cvScalar(0,0, 128) ]

OUT = cStringIO.StringIO()

class Edge:
    def __init__(self):
        self.start = None 
        self.end = None 
        self.points = []
    def addPoint(self, x):
        self.points.append(x)

class Linker:

    #find the corners of images, and save all corner points in self.vKeyPoints
    def __findedge(self, filename): 
        tmpimg = highgui.cvLoadImage (filename)
        self.img = cv.cvCreateImage(cv.cvSize(int(tmpimg.width * self.enlarge), int(tmpimg.height * self.enlarge)), 8, 3)
        cv.cvResize(tmpimg, self.img, cv.CV_INTER_LINEAR)
        if (self.drawimage):
            self.drawimg = cv.cvCloneImage(self.img)
        else:
            self.drawimg = cv.cvCreateImage(cv.cvGetSize(self.img), 8, 3)
        greyimg = cv.cvCreateImage(cv.cvSize(self.img.width, self.img.height), 8,1)
        cv.cvCvtColor(self.img, greyimg, cv.CV_BGR2GRAY)
        self.allcurve = []
        for i in range(80, 200, 20):
            bimg = cv.cvCloneImage(greyimg) 
            cv.cvSmooth(bimg, bimg, cv.CV_MEDIAN, 9)
#            cv.cvSmooth(bimg, bimg, cv.CV_BILATERAL, 9)
#            cv.cvSmooth(bimg, bimg, cv.CV_BLUR, 9)
#            cv.cvSmooth(bimg, bimg, cv.CV_BLUR, 9)
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
        if (self.drawcontour):
            cv.cvDrawContours (self.drawimg, cont, _white, _white, 1, 1, cv.CV_AA, cv.cvPoint (0, 0))
        idx = 0
        for c in cont.hrange():
            PointArray = cv.cvCreateMat(1, c.total, cv.CV_32SC2)
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

    def __edgededup(self):
        matcher = NW.NWMatch()
        todelete = {}
        for i, e1 in enumerate(self.edges):
            if (len(e1.points) < 5):
                todelete[i] = 1
                continue
            for j, e2 in enumerate(self.edges):
                 if (j <= i ) : continue
                 if (todelete.has_key(i)): continue
                 if (len(e2.points) < 15): 
                      todelete[j] = 1
                      continue
                 if (e2.start == e1.start and e2.end ==  e2.end):
                      x = [[e.x, e.y] for e in e1.points]
                      y = [[e.x, e.y] for e in e2.points]
                      #print x, y
                      cost = matcher.Align(x, y)
                      if (cost < 2 * max(len(x), len(y))) : 
                          todelete[j] = 1
        print "number of edge before dedup %d " % ( len(self.edges)) 
        self.edges = [e for i, e in enumerate(self.edges) if not todelete.has_key(i)]  
        print "number of edge after dedup %d " % ( len(self.edges)) 
                      
         
    def __evenSample(self, npoint):
        totalpoints = 0
        for e in self.edges:
            totalpoints += len(e.points)
        ndist = totalpoints / (0.001 + npoint)
        self.allselected = []
        for c in self.edges:
            c.selected = []
            c.selected.append(self.points[c.start])
            for pi, p in enumerate(c.points):
                if (int((pi + 1) / ndist) - int(pi / ndist)) == 1:
                    c.selected.append(cv.cvPoint(int(p.x), int(p.y)))
            c.selected.append(self.points[c.end])

    def __link(self):
        myfont = cv.cvInitFont(cv.CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5)
        kkk = 0
        self.edges = []
        for curve in self.allcurve:
            showpt = []
            state = 0  
            currentPoint = None 
            currentPointIdx = -1
            cumulate = 0
            dcurve = curve + curve
            curlen = len(curve)
            ptcount = 0
            pointseq = []
            for c in dcurve:
                if (ptcount > curlen): break
                cumulate += 1
                for kk in range(len(self.points)):
                     k = self.points[kk]
                     if (abs(c.x - k.x) + abs(c.y - k.y) < 15):
                          if (currentPoint != k or cumulate > 40):
                              state += 1
                              currentPoint = k
                              currentPointIdx = kk
                              cumulate = 0
                              pointseq.append(kk)
                if (state > 0):
                     showpt.append([c, state, currentPointIdx])
                     ptcount += 1
            if (state > 1):
                kkk += 1
                cnt = 0
                pret = -1
                e = None 
                for s,t, cp in showpt:
                    if (cp != pret):
                        if e != None:
                             e.end = cp
                        e = Edge()
                        self.edges.append(e)
                        e.start = cp 
                        pret = cp 
                    cnt += 1
                    if (t < state):
                         e.addPoint(s)
                         #print  "%d\t%3.2f\t%3.2f\t%d\t%d\t%d" % (kkk, s.x, s.y, cnt, pointseq[t - 1], pointseq[t])
                e.end = showpt[-1][2]
        print >> OUT, "seq\tptn\tx\ty\t" 
#        self.__edgededup()
        self.__evenSample(self.npoints)
        for ie, e in enumerate(self.edges):
            print  "P(%d) <-> P(%d) length %d, selected %d" % (e.start, e.end, len(e.points), len(e.selected))
            for d in e.points:
                cv.cvSet2D(self.drawimg, int(d.y), int(d.x), color[3])
            for ip, p in enumerate(e.selected):
                cv.cvDrawCircle(self.drawimg, p, 2, cv.cvScalar(255,255,0,0))
                print >> OUT, "%d\t%d\t%d\t%d" % (ie, ip, p.x, p.y) 
               
#            for id in range(1, len(e.selected)):
#                cv.cvLine(self.drawimg, e.selected[id-1], e.selected[id], cv.CV_RGB(255,0,0))
            #highgui.cvShowImage ("Corner1", self.drawimg)
            #highgui.cvWaitKey (0)

                             
        
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


    def LinkPoints(self, fname):
        self.__findedge(fname)
        self.DrawKeyPoints()
        self.__link()

    def HarrisPoints(self, imgfile):
        self.points = []
        self.drawimg = highgui.cvLoadImage (imgfile)
        c = 1
        try:
            gray = cv.cvCreateImage (cv.cvGetSize (self.drawimg), 8, 1)
            cv.cvCvtColor(self.drawimg, gray, cv.CV_BGR2GRAY)
            eig = cv.cvCreateImage (cv.cvGetSize (self.drawimg), 32, 1)
            tmpimg = cv.cvCreateImage (cv.cvGetSize (self.drawimg), 32, 1)
            p =   cv.cvGoodFeaturesToTrack(gray, eig, tmpimg, 100, 0.1, 20, None, 7, 1, 0.04 )
            for x in p:
                 cv.cvCircle( self.drawimg, x, 3, cv.CV_RGB(0,255,0), 8, 0 );
                 self.points.append(x)

        except Exception,e:
            print e
            print 'ERROR: problem handling '+ imgfile 
 
    
    def LoadPoints(self, ptfile):
        self.points = []
        for line in fileinput.input(ptfile):
            dr = line.strip("\n").strip("\r").split(" ")
            ds = [d.strip(" ") for d in dr]
            dt = [d for d in dr if d != ""]
            x = float(dt[0]) * self.enlarge
            y = float(dt[1]) * self.enlarge 
            pt = cv.cvPoint(int(x), int(y))
            self.points.append(pt)
              
    def __init__(self, drawcontour, drawimage, enlarge, npoints):
        self.output  = None
        self.bDraw  = 0
        self.bFilter = 1
        self.drawimg = None
        self.drawcontour = drawcontour
        self.drawimage = drawimage
        self.enlarge = enlarge
        self.npoints = npoints
        
    
    
def main():
    usage = "%s [options]  <imgfile> " % (sys.argv[0])
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--display', action="store_true", dest = 'display', default = False, help = 'display the image')
    oparser.add_option('-c','--contour', action="store_true", dest = 'contour', default = False, help = 'show object contour')
    oparser.add_option('-i','--image', action="store_true", dest = 'image', default = False, help = 'show original images')
    oparser.add_option('-n', '--number', dest = 'num', type='int', default = 200 , help = 'the number of feature points')
    oparser.add_option('-x','--enlarge', dest = 'enlarge', default = 1.0 , type = float,  help = 'resize images, default:1.0')
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')
    oparser.add_option('-p', '--pointfile', dest = 'pointfile', default = None, help = 'use pointfile ')
    oparser.add_option('-r', '--harris', dest = 'harris', default = False, action = "store_true", help = 'use harris detector')
    oparser.add_option('-s', '--save', dest = 'save', default = None, help = 'save the img file')
    

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 2:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    if (options.pointfile == None and options.harris == None): 
        print >> sys.stderr, "either of  pointfile and harris can be valid"
        sys.exit(1)

    highgui.cvNamedWindow ("Corner1", 1)
    ct = Linker(options.contour, options.image, options.enlarge, options.num)
    if (options.pointfile): 
        ct.LoadPoints(options.pointfile)
        ct.LinkPoints(args[1])
    else:
        ct.HarrisPoints(args[1])
        ct.LinkPoints(args[1])
    highgui.cvShowImage ("Corner1", ct.drawimg)
    highgui.cvWaitKey (0)   
    if (options.save):
        highgui.cvSaveImage(options.save, ct.drawimg)    
    if (options.output):
        f = open(options.output, "w")
        f.write(OUT.getvalue())
        f.close()
        OUT.close()

if __name__ == '__main__':
    # load the image gived on the command line
    main()
