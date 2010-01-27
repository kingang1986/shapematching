#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct, pickle
from makeIdx import MSIndexer


def main():
  
    usage = "%prog [options] <datafile.idxed> <datafile.query>"
    version = "%prog 0.2\n Query an Indexed MSS database. \nLongbin Chen, longbinc@yahoo.com"

    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-k', '--keywordlength', dest = 'keywordlength', type='int',default = 8 , help = 'keyword length, default: 8')
    oparser.add_option('-c', '--codebook', dest = 'codebook', default = None , help = 'codebook file,  default: None')
    oparser.add_option('-d', '--database', dest = 'database', default = None, help = 'use database')
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    idxer = MSIndexer()
    idxer.loadIdx(args[1])
    idxer.keywordlength = options.keywordlength

    # search
    fin = open(args[2], "r")
    while 1:
            line = fin.readline()
            if not line: break;
            line = line.strip("\n\r\t ")
            if (len(line) == 0): continue
            if (line[0]  == ">"):
                cmm = line 
            else:
                seq = line.split(" ")
                print >> sys.stderr, "search for  %s " % cmm
                idxer.queryDB(seq)
    fin.close()
       
       
        
    
if __name__ == '__main__':
    main()
