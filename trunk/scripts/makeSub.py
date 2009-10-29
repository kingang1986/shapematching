#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct, array
import SW

def GetSubScore( index1, index2):
    f = 0.0
    for k in range(len(index1)):
        d  =  index1[k] - index2[k]
        m  = abs(index1[k]) + abs(index2[k])
        f = f + abs(d)
    f = 1024 - f * 200 
    if (f < 0):
         return -6
    f =  math.log(f) / math.log(2) - 6
    return f

def GetSubScore_bak( index1, index2):
    f = 0.0
    for k in range(len(index1)):
        d  =  abs(index1[k] - index2[k])
        m  =  abs(index1[k]) + abs(index2[k])
        if ( m > 0):
            f += d  / m
    f = 5 - f
    return f

def fifp(value):
    str = ""
    value = int(value)
    if (value == 5):
        str = "  %d" % (5) 
    else:
        str = " %d" % (-1)
    return str
def ffp(value):
    str = ""
    value = int(value)
    if (value < 0):
        str = " %d" % (value) 
    else:
        str = "  %d" % (value)
    return str
def main():
  
    usage = "%prog [options]  <codebook.cbk> <sub_matrix_file>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)

    (options, args) = oparser.parse_args(sys.argv)

    options.unknown = -6
    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    fileobj = open(args[1], 'rb')
    magicnum = fileobj.readline()
    line2 = fileobj.readline()
    N, dim = [int(x) for x in line2.strip(" ").split()]
    print magicnum,
    print "%d lines, each line %d dim "% (N, dim) 
    res = []
    for i in range(N):
        v = []
        for j in range(dim):
            t = struct.unpack('>f', fileobj.read(4))
            v.append(t[0])
        res.append(v)
    fileobj.close()

    fout = open(args[2], "w")
    for i in range(N):
        for j in range(N):
             sub = GetSubScore(res[i], res[j]) 
             t = int(sub* 25)
             fout.write("%s," % (ffp(sub)))
        fout.write(" -1,%s,\n" % ffp(options.unknown)) 
    for i in range(N + 1):
        fout.write(" -1,") 
    fout.write("%s,\n" % ffp(options.unknown))
    for i in range(N+1):
        fout.write("%s," % ffp(options.unknown)) 
    fout.write("  1\n")

    fout.close()
    
if __name__ == '__main__':
    main()
