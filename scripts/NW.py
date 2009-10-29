#!/usr/bin/env python
"""NWMatch - a Python class for implementation of Needleman-Wunsch algorithm
"""

__author__ = "Longbin Chen(longbinc@yahoo.com)"
__version__ = "$Revision: 1.0 $"
__history__="""
    20090831:  Created

    """

import sys, os, glob,  optparse, fileinput, string, pickle, Numeric, math
from MSS import * 


OP_MATCH = 1
OP_SUBST = 2
OP_DELET = 3 
OP_INSRT = 4 
NO_OP =  -1
NO_OP_BOTH = -1
NO_OP_A =   -2
NO_OP_B = -3 

class NWMatch:

    def SetModel(self, weights, matchmap, gapmap):
        self.weights = weights[:]
        self.matchmap = matchmap[:]
        self.gapmap = gapmap[:]
        
    def GetSubCost(self, index1, index2):
    	f = 0.0
        for k in self.keys:
    	    d  = (self.seq1.points[index1].getvalue(k) -  self.seq2.points[index2][k])  
    	    m  = (self.seq1.points[index1][k] +  self.seq2.points[index2][k])  
            if (m > 0):
    	        f += d * d / m 
        return  f 
 
    def Align(self, seq1, seq2, delcost = 2):
        self.delCost = delcost
    	self.seq1 = seq1
    	self.seq2 = seq2
        self.keys = [k for k in self.seq1.points[0].getvalues().keys() if k.startswith("f")]
    	MatchX = []
    	MatchY = []
        matchscores = []
        m = len(seq1.points) 
        n = len(seq2.points)
        cost = Numeric.zeros((m + 1, n + 1), 'f')
        path = Numeric.zeros((m + 1, n + 1), 'i')
        for i in range(m + 1):
            cost[i][0] = self.delCost * i
            path[i][0] = OP_DELET
        for i in range(n + 1):
            cost[0][i] = self.delCost * i
            path[0][i] = OP_INSRT
        path[0][0] = 0
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                i1 = i - 1
                j1 = j - 1
                fSubstitute = cost[i1][j1] + self.GetSubCost(i1, j1) 
                fDelete = cost[i1][j] + self.delCost
                fInsert = cost[i][j1] + self.delCost
                if (fSubstitute < fDelete and fSubstitute < fInsert):
                    cost[i][j] = fSubstitute
                    path[i][j] = OP_SUBST 
                elif (fSubstitute > fDelete and fDelete < fInsert):
                    cost[i][j] = fDelete
                    path[i][j] = OP_DELET
                else:
                    cost[i][j] = fInsert 
                    path[i][j] = OP_INSRT
        aligncost = cost[m][n]
    #    for x in range(m + 1):
    #        for y in range(n + 1):
    #            print path[x][y],
    #        print 
        u = m
        v = n
        while u > 0 or v > 0:
            if path[u][v] == OP_SUBST:
               MatchX.append(u - 1)
               MatchY.append(v - 1)
               u -= 1
               v -= 1
            elif path[u][v] == OP_DELET:
               MatchX.append(u - 1)
               MatchY.append(- 1)
               u -= 1
            else: 
               MatchX.append(- 1)
               MatchY.append(v - 1)
               v -= 1
        return aligncost, MatchX, MatchY



def main():

    usage = "%prog [options] <feature_file1> <feature_file2>\nMatch two sequence using Needleman Wusche algorithm"
    version = "%prog 0.2\nLongbin Chen, 08/31/2009,  longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-o', '--output', dest = 'output', default = None, help = 'output file')
    oparser.add_option('-m', '--model', dest = 'model', default = None, help = 'model file') 


    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    mss1 = MSS()
    mss1.load(args[1])
    mss2 = MSS()
    mss2.load(args[2])

    ######## Typical usage of Needleman-Wunsch algorithm
    ######## seq1, seq2 are two lists 
    matcher = NWMatch()
    cost, matchX, matchY = matcher.Align(mss1.seqs[0], mss2.seqs[0])
    print cost
    if (options.output):
       f = open(options.output, "w")
       for i in range(len(matchX)):
           f.write("%d\t%d\n" %(matchX[i], matchY[i]))
       f.close()
    ######## End of Typical usage of Smith Waterman algorithm

     	
if __name__ == '__main__':
    main()
