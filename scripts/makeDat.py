#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct



    
def main():
  
    usage = "%prog [options] <datafile.tsv> <datafile.dat>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"

    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-e', '--exclude-feature', dest = 'exclude', default = None, help = 'exclude features from .tsv file')
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    # input header
    ftsv  = open(args[1], "r")
    header = ftsv.readline().strip("\n\t").split("\t") 
    featidx = range(len(header)) 
    if (options.exclude):
        featidx = [idx for idx in featidx if header[idx] not  in options.exclude.split(",")] 
    count = 0
    fmax = -1e10
    fmin = 1e10
    while 1:
         line = ftsv.readline()
         if (not line): break;
         t = [float(a) for a in line.strip("\n\t").split("\t")]
         mm = max(t)
         if (fmax < mm): fmax = mm
         mn = min(t)
         if (fmin > mn): fmin = mn
         count += 1
    ftsv.close()
    ftsv = open(args[1], "r")
    ftsv.readline()  # header 
    
    fout = open(args[2], "wb")
    fout.write("DAT0.58\n")
    fout.write("%d %d\n"% (count, len(featidx)))
    fout.write("%f %f\n"% (fmin, fmax))
    while 1:
        line = ftsv.readline()
        if (not line): break;
        linedat = [float(t) for t in line.strip('\n\t').split('\t')]
        for i in featidx:
            fout.write(struct.pack('>f', linedat[i]))
    ftsv.close()
    fout.close()
        
    
if __name__ == '__main__':
    main()
