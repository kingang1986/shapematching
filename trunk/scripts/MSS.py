#!/usr/bin/env python
'''
 The class to read, write the multi-sequence shape
 ------------
 Created by longbinc@yahoo.com
'''

import sys, math, string, optparse, fileinput, cStringIO, random


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

class Sequence:
    def __init__(self):
        self.points = []
        
    def addPoint(self, newpt):
        self.points.append(newpt)
    def __getitem__(self, k):
        return self.points[k]
        
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
        
        feaidx = [header.index(fld) for fld in header if fld.startswith("f")]
        self.featname = [fld for fld in header if fld.startswith("f")]
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

    def save(self, filename):
        res = self.seqs[0].points[0].getvalues().keys()
        res.sort()
        self.featname = [f for f in res if f.startswith("f") and f not in ['seq', 'pnt', 'x', 'y']]
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

