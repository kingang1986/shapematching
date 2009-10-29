#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct, pickle
from makeIdx import MSIndexer


def main():
  
    usage = "%prog [options] <datafile.idxed> <datafile.query>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"

    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-d', '--database', dest = 'database', default = None, help = 'use database')
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    idxer = MSIndexer()
    idxer.loadIdx(args[1])

    # search
    fin = open(args[1], "r")
    while 1:
        line = fin.readline()
        if not line: break;
        line = line.strip("\n\r\t ")
        if (len(line) == 0): continue
        if (line[0]  == ">"):
            cmm = line 
        else:
            print "search result for %s " % (cmm)
            seq = line.split(" ")
            idxer.findSequence(seq)
    fin.close()
    idxer.saveIdx(args[2])
    
       
        
    
if __name__ == '__main__':
    main()
