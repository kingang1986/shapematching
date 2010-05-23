#!/usr/bin/env python
'''
 The class to read, write the multi-sequence shape
 ------------
 Created by longbinc@yahoo.com
'''

import sys, math, string, optparse, fileinput, cStringIO, random
from opencv import cv
from opencv import highgui


class myPoint:
    def __init__(self):
        self.x = 0
        self.y = 0
        self.index = -1
        
    def getvalues(self):
     	 return self.__dict__

    def getvalue(self, key):
        return self.__dict__[key]
    def __getitem__(self, k):
        return self.__dict__[k]

    def addvalue(self, name, value):
         self.__dict__[name] = value
    def getCvPoint(self):
        return cv.cvPoint(int(self.x), int(self.y))

class Sequence:
    def __init__(self):
        self.points = []
        
    def addPoint(self, newpt):
        self.points.append(newpt)

    def __getitem__(self, k):
        return self.points[k]

    def paint(self, img):
        for p in self.points:
            cv.cvDrawCircle(img, p.getCvPoint(), 2, cv.cvScalar(0, 0, 255,0))
        for i in range(len(self.points) - 1):
            cv.cvLine(img, self.points[i].getCvPoint(), self.points[i + 1].getCvPoint(), cv.cvScalar(255,255,255,0), 1)
    def cut_seq(self, pos1, pos2): # remove the sequence from pos1 to pos2, and return two (rest) sequence 
        s1 = Sequence()
        s2 = Sequence()
        for i in range(pos1):
            s1.points.append(self.points[i])
        for j in range(pos2,len(self.points)):
            s2.points.append(self.points[j])
        return s1, s2
        
class MSS:
    def __init__(self):
        self.seqnum = 0
        self.seqs = []
        self.featname = []
        
    def load(self, filename):
        self.seqs = []
        idx = 0
        fp = open(filename, "r") 
        header = fp.readline().strip("\n").split("\t")
        
        feaidx = [header.index(fld) for fld in header if fld.startswith("f") or fld.startswith("g")]
        self.featname = [fld for fld in header if fld.startswith("f") or fld.startswith("g")]
        seq_idx = header.index("seq")
        pnt_idx = header.index("ptn")
        x_idx = header.index("x")
        y_idx = header.index("y")
        seq = Sequence()
        pre_seq = -1
        for line in fp.readlines():
            pt = myPoint()
            line = line.strip("\n")
            dr = line.split('\t')
            pt.x = string.atof(dr[x_idx])
            pt.y = string.atof(dr[y_idx])
            for f in feaidx: 
                pt.addvalue(header[f], string.atof(dr[f]))
            if (dr[seq_idx] != pre_seq and pre_seq != -1):
                self.seqs.append(seq)
                seq = Sequence()
            pre_seq = dr[seq_idx]
            seq.addPoint(pt)
        self.seqs.append(seq)

    def feature_filter(self, featurenames):
        res = []
        for f in featurenames:
            if (not(f in ['seq','pnt','x','y'])):
                res.append(f)
        return res

    def paint(self, img): # draw MSS on img
        for s in self.seqs:
            s.paint(img)
        

    def save(self, filename):
        if len(self.seqs) == 0: return
        if len(self.seqs[0].points) == 0: return
        res = self.seqs[0].points[0].getvalues().keys()
        res.sort()
        self.featname = [f for f in res if (f.startswith("f") or f.startswith("g")) and f not in ['seq', 'pnt', 'x', 'y']]
        fout = open(filename, "w")
        fout.write("seq\tptn\tx\ty\t" + "\t".join(self.featname) + "\n")
        for si, seq in enumerate(self.seqs):
            for pi, pnt in enumerate(seq.points):
                var = pnt.getvalues()
                fout.write("%d\t%d\t%3.2f\t%3.2f\t%s\n" % (si, pi, pnt.x, pnt.y, "\t".join(["%3.4f" % (var[f]) for f in self.featname])))
        fout.close()    
    
    def getAllPoints(self):
        allpt = []
        for seq in self.seqs:
            for pnt in seq.points:
                allpt.append(pnt)
        return allpt

def main():
    usage = "%s [options]  <shape file>" % (sys.argv[0])
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 2:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    mss = MSS()
    mss.load(args[1])
    mss.save(args[1] + ".cpy")

    
if __name__ == '__main__':
    # load the image gived on the command line
    main()

