#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct, array

    
def main():
  
    usage = "%prog [options]  <datafile> <mssfile>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-t', '--template', dest = 'template', default = None, help = 'template file ')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    fileobj = open(args[1], 'rb')
    fout = open(args[2], "w")
    magicnum = fileobj.readline()
    line2 = fileobj.readline()
    N, dim = [int(x) for x in line2.strip(" ").split()]
    line3 = fileobj.readline()
    min, max= [float(x) for x in line3.strip(" ").split()]
    if (options.template):
       ftempl = open(options.template, "r")
       header = ftempl.readline()
       fout.write(header)
    print magicnum,
    print "%d lines, each line %d dim "% (N, dim) 
    print "%f min, %f max "% (min, max) 
    for i in range(N):
        if (options.template):
            dr = ftempl.readline().strip("\n\t").split("\t")
            fout.write("\t".join(dr[0:4]))
        for j in range(dim):
            t = struct.unpack('>f', fileobj.read(4))
            fout.write("\t%3.5f" % t[0])
        fout.write("\n")
    if (options.template):
        ftempl.close()
    fileobj.close()
    fout.close()
    
if __name__ == '__main__':
    main()
