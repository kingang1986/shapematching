#!/bin/python 
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
from scipy import * 
import math
import getopt
import fileinput
import string   
import pickle

        
OP_MATCH = 1
OP_SUBST = 2
OP_DELET = 3 
OP_INSRT = 4 
NO_OP =  -1
NO_OP_BOTH = -1
NO_OP_A =   -2
NO_OP_B = -3 
#w = [0.006563, 0.004539, 0.000521, 0.000008, 0.005244, 0.009312, 0.005698, 0.010103, 0.017570, 0.005054, 0.016251, 0.017254, 0.000000,  0.013414, 0.012713]
#w = [-0.018364, 0.018848, 0.041341, -0.000004, 0.026212, 0.035691, 0.028115, 0.041407, 0.170188, 0.041351, 0.111188, 0.126718, 0.000000,  0.042128, 0.129607]
#w = [0.036445, -0.052975, -0.064812, 0.059437, -0.041999, -0.056073, -0.042927, -0.052255, -0.171573, -0.046712, -0.101602, -0.122834, 0.000000,  -0.069309, 0.109205]
#w  =[-0.041190, -0.080295, -0.134056, -0.022961, 0.724615, -0.220173, -0.472291, -0.122831, 0.014997, -0.133782, -0.212366, -0.270956, 0.000000,  -0.080177, 0.231213]
#w  =  [-0.013253,-0.066678,-0.085552,-0.023620,0.588756,-0.221439,-0.298726,-0.128879,-0.018632,-0.117025,-0.098258,-0.137467,0.000000,-0.033016,0.165622]

class SmithWaterman:

    def LoadModel(self, modelname):
    	fmodel = fopen(modelname, 'r')
    	line1 = fmodel.readline()
    	line2 = fmodel.readline()
    	if (line2[-1] == '\n'):
    		line2 = line2[:-1]
    	dr = line2.split(' ')
    	self.w = [ string.atof(x) for x in dr] 
    	self.feature_length  = len(dr)
    	
    	
    def InitModel(self, feature_length):
        print "initialize model %d"%(feature_length)
    	self.w = [0] * ( 2 * feature_length + 2)
    	for i in range(feature_length):
    	    self.w[i] = -0.1
    	self.w[4] = -0.1
    	self.w[feature_length] = 0.2
    	self.w[feature_length * 2 + 1] = -0.2
    	self.feature_length  = feature_length

    def GetSubScore(self, index1, index2):
    	f = 0.0
    	for i in range(self.feature_length):
    	    d  = (self.seq1[index1][i] -  self.seq2[index2][i])  
    	    m  = (self.seq1[index1][i] +  self.seq2[index2][i])  
    	    if (m > 0):
    	    	f += self.w[i] * d * d / m
        return  self.w[self.feature_length] + f 

    def GetGapScore(self, indexi):
    	for i in range(self.feature_length):
    	#    f +=  self.w[i + self.feature_length + 1] * abs(self.seq1[indexi][i]) 
        #return self.w[self.feature_length * 2 + 1] + f

    def GetInsertCost(self, indexi, indexj):
    	return  0.1
    	#for i in range(self.feature_length):
    	#    f +=  self.w[i + self.feature_length + 1] * abs(self.seq2[indexj][i]) 
        #return self.w[self.feature_length * 2 + 1] + f
 
    def PrintPath(self, m):
        for L in m:
           for k in L:
               print "%3.3f "%(k),
           print ""              
     
    def Align(self, seq1, seq2):
    	self.seq1 = seq1
    	self.seq2 = seq2
    	self.feature_length = len(seq1[0])
    	MatchX = []
    	MatchY = []
        m = len(seq1) 
        n = len(seq2)
        cost = zeros((m + 1, n + 1), float)
        path = zeros((m + 1, n + 1), int)
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                i1 = i -1
                j1 = j -1
                fSubstitute = cost[i1][j1] + self.GetSubScore(i1, j1) 
                fDelete = cost[i1][j] + self.GetDeleteCost(i1, j)
                fInsert = cost[i][j1] + self.GetInsertCost(i, j1)
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
        npos = argmax(cost)
        nRow = npos / (n + 1)
        nCol = npos % (n + 1)
        maxrow =nRow
        maxcol = nCol
        alignMatch = []
        subsvalue = []
        while (nRow > 0 and nCol > 0 ):
            t = path[nRow][nCol]
            if (t == 1):
                MatchX.append(nRow - 1)
                MatchY.append(nCol - 1)
                nRow -= 1
                nCol -= 1
                alignMatch.append(OP_SUBST)
            elif (t == 2):
                nRow -= 1
                alignMatch.append(OP_DELET)
            elif (t == 3):
                nCol -= 1
                alignMatch.append(OP_INSRT)
            elif (t == 0):
                break
        align = []
        if (nRow > nCol):
            for i in range(nRow - nCol):	
                align.append(NO_OP_A)
        else:
            for i in range(nCol - nRow):	
                align.append(NO_OP_B)
        for i in range(min(nCol, nRow)):
        	align.append(NO_OP)
        alignMatch.reverse()
        align = align[:] + alignMatch[:]
        return MatchX, MatchY
        
 
        
def loaddata(filename, headerline):
    res = []
    line_count = 0
    for line in fileinput.input(filename):
    	if (line[-1] == '\n'):
    		line = line[:-1]
        if (len(line) < 2):
        	continue
        if (line_count == 0) and (headerline == 1):
       	    line_count += 1
            continue
        else:
       	    line_count += 1
            dr = line.split('\t')
            line_res = []
            for d in dr:
                line_res.append(string.atof(d))
            res.append(line_res)
    return res
 

def usage():
    print "usage: %s [-m modelfile] [-o outputfile] [-h --header] datafile1 datafile2"%sys.argv[0]
    
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:o:h", ["modelfile=", "output=", "--header"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    modelfile = None     
    output = None
    headerline = 1
    for o, a in opts:
        if o in ("-o", "--output"):
            output = a
        if o in ("-m", "--modelfile"):
            modelfile = a
        if o in ("-h", "--header"):
            headerline = 0
    if (len(args)) != 2:
        usage()
        sys.exit(2)
    seq1 = loaddata(args[0], headerline)
    seq2 = loaddata(args[1], headerline)
    matcher = SmithWaterman()
    if (modelfile):
    	matcher.LoadModel(modelfile)
    else:
        matcher.InitModel(len(seq1[0]))
    cost,align = matcher.Align(seq1, seq2)
    if (output):
	foutput = open(output,"w")
	for k in align[:-1]:
	    foutput.write("%d,"%(k))
	foutput.write("%d"%(align[-1]))
	foutput.close()

     	
if __name__ == '__main__':
    # load the image gived on the command line
    main()
