#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct


    
def main():
  
    usage = "%prog [options] <ilist> "
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"

    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-u', '--uniq', dest = 'uniq', action = "store_true", default = True, help = 'uniq pairs, default False')
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 2:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    f = open(args[1], "r")
    existline = []
    while 1:
       line = f.readline() 
       if (not line): break;
       line = line.strip("\n")
       for i, e in enumerate(existline):
           print e +"\t"+line+("\t%d\t%d" %(i, len(existline)))
           if (not options.uniq):
               print line + "\t" + e ("\t%d\t%d" %(len(existline), i ))
       existline.append(line)
    f.close()
       
        
        
    
if __name__ == '__main__':
    main()
