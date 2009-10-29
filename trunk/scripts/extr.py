#!/usr/bin/env python
__author__ = "Longbin Chen  lbchen@cs.ucsb.edu"
__version__ = "$Revision: 1.0 $"
__history__="""
"""

import sys, os, optparse, fileinput, math, string 
from CntPoint import *
from CntDAS import *
from CntAngle import *
from CurvePoint import * 
from CntSC import *

def filter_img(fname):
    res = []
    for f in fname:
        if ( lower(f[-4:]) == '.jpg') or ( lower(f[-4:]) == '.bmp'):
            res.append(f)
    return res

def extract_dir(inputdir, outputdir, options):
    if (inputdir[-1] == "/"):
       inputdir = inputdir[:-1]
    if (outputdir[-1] == "/"):
       outputdir = outputdir[:-1]
    ct = CurvePoint()
    sc = CntSC()
    das = CntDAS()
    agl  = CntAngle()
    allfile = os.listdir(inputdir)
    to_file = options.ter
    start_file = options.start
    if (to_file == -1):
        to_file = len(allfile)
    if (to_file > len(allfile)):
        to_file = len(allfile)
    if (start_file < 0):
       start_file = 0

    for i in range(start_file, to_file):
        f = allfile[i]
        inputfile = inputdir + "/"+ f
        outputfile = outputdir + "/" + f
        ct.LoadCont(inputfile)
        for c in ct.allselected:
            if (options.shapecontext):
                sc.ExtractFeature(c)
            if (options.angle):
                agl.ExtractFeature(c, 0)
            if (options.das):
                das.ExtractFeature(c, 0)
        ct.Save(outputfile)

        

def extract_file(input, output, options):
    ct = CurvePoint()
    sc = CntSC()
    das = CntDAS()
    agl  = CntAngle()
    
    ct.LoadCont(input)
    for c in ct.allselected:
        if (options.shapecontext):
           sc.ExtractFeature(c)
        if (options.angle):
           agl.ExtractFeature(c, 0)
        if (options.das):
           das.ExtractFeature(c, 0)
    ct.Save(output)


        
if __name__=='__main__':
 
    usage = "%prog [options] <img/img_dir> <featurefile/feature_dir>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--directory', action="store_true", dest = 'directory', default = False, help = 'input and output are directories, extract feature for each image in the input directory and save the feature files in the output directory')
    oparser.add_option('-n', '--number', dest = 'num', type='int',default = 200 , help = 'the number of feature points')
    oparser.add_option('-s', '--start', dest = 'start', type='int',default = 0 , help = 'start from image #, default 0, from the first one, valid only when -d')
    oparser.add_option('-t', '--terminate', dest = 'ter', type='int',default = -1 , help = 'terminate with image #, default to the end')
    oparser.add_option('-e', '--even', action="store_false",dest = 'even', default = True,  help = 'evenly sample feature points from the contour(default true)')
    oparser.add_option('-a', '--angle', action="store_false",dest = 'angle', default = True,  help = 'extract angle features(default true)')
    oparser.add_option('-c', '--shapecontext', action="store_false",dest = 'shapecontext', default = True,  help = 'extract shape context  features(default true)')
    oparser.add_option('-w', '--das', action="store_false",dest = 'das', default = True,  help = 'extract DAS features(default true)')

    (options, args) = oparser.parse_args(sys.argv)
    
    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
 
    if (options.directory):
       extract_dir(args[1], args[2], options) 
    else:
       extract_file(args[1], args[2], options)
 
            
            
