#!/usr/bin/python
"""mpeg7_feature.py - a Python class
"""

__author__ = "Longbin Chen  lbchen@cs.ucsb.edu"
__version__ = "$Revision: 1.0 $"
__history__="""

"""
import sys
import getopt
import os
from string import *
import fileinput

headline = "idsc 0	idsc 1	idsc 2	idsc 3	idsc 4	idsc 5	idsc 6	idsc 7	idsc 8	idsc 9	idsc10	idsc11	idsc12	idsc13	idsc14	idsc15	idsc16	idsc17	idsc18	idsc19	idsc20	idsc21	idsc22	idsc23	idsc24	idsc25	idsc26	idsc27	idsc28	idsc29	idsc30	idsc31	idsc32	idsc33	idsc34	idsc35	idsc36	idsc37	idsc38	idsc39	idsc40	idsc41	idsc42	idsc43	idsc44	idsc45	idsc46	idsc47	idsc48	idsc49	idsc50	idsc51	idsc52	idsc53	idsc54	idsc55	idsc56	idsc57	idsc58	idsc59"

def filter_img(fname):
    res = []
    for f in fname:
        if ( lower(f[-4:]) == '.jpg') or ( lower(f[-4:]) == '.bmp'):
            res.append(f)
    return res

def filter_txt(fname):
    res = []
    for f in fname:
        if ( lower(f[-4:]) == '.txt') :
            res.append(f)
    return res

def modify(fname):
    fp = open(fname, 'r') 
    lines = fp.readlines() 
    fp.close()
    fq = open(fname, 'w')
    fq.write(headline)
    fq.write('\n')
    for lin in lines:
        fq.write(lin)
    fq.close()
    



def usage():
    print "usage: %s [-h] [-s startfrom] [-t to] filefolder  "%sys.argv[0]

if __name__=='__main__':
 
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:t:", ["help", "startfrom", "to"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    if (len(args)) != 1:
        usage()
        sys.exit(2)      

    npoint = 20
    start = 0
    to = -1
     
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-s", "--startfrom"):
            start = string.atoi(a)
        if o in ("-t", "--to"):
            to = string.atoi(a)
    allfile = filter_txt(os.listdir(args[0]))
    if (to == -1):
        to = len(allfile)
    for i in range(start, to):
    	f = allfile[i]
        g = args[0] + "\\"+ f
        modify(g)
  
 
            
            