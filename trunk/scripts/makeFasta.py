#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct

def main():
  
    usage = "%prog [options]  <datafile.chn> <datafile.tsv> <datafile.fasta>"
    version = "%prog 0.2\nGenerate fasta format data using codebook and original TSV feature file(to provide sequence info). Longbin Chen, longbinc@yahoo.com"

    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-l', '--label', dest = 'label', default = None, help = 'label file: default: None')
    oparser.add_option('-C', '--char', dest = 'char', action="store_true", default = False, help = 'Represent cluster with characters. default: False')
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 4:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    # input header
    outputstring = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']
    strlabel = ""
    
    if (options.label):
        idx = int(args[1].split("/")[-1])
        flabel = open(options.label, "r")
        r = ""
        for i in range(idx):
            r = flabel.readline().strip("\r\n") 
        strlabel = "| %s" % (r)
        flabel.close()
    fout = open(args[3], "w")
    seqn = 0
    fout.write(">mn| %s | %d %s\n" % (args[1], seqn, strlabel))
    fdat  = open(args[1], "r")
    fall = open(args[2], "r")
    fall.readline()
    fdat.readline() # magic number line
    num = int(fdat.readline().strip("\n\t"))
    fdat.readline() # code book size
    for i in range(num):
       seq = int(fall.readline().split("\t")[0])
       if seq != seqn:
           fout.write("\n>mn| %s | %d %s\n" % (args[1], seq, strlabel))
           seqn = seq
       d = int(fdat.readline().strip("\n\t"))
       if (options.char):
           fout.write(outputstring[d-1])
       else:
           fout.write("%d " % d)
    fout.write("\n")
    fout.close()
    fdat.close()
    
if __name__ == '__main__':
    main()
