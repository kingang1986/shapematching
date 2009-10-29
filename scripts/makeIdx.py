#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct, pickle



class  MSIndexer:
    def __init__(self):
        self.keydic = {}
        self.seqdic = {}
        self.seqinfo = [] 
        self.seqdata = [] 
        self.windowsize = 21 
        self.seqnumber = 0

    def __keyDistance(self, k1, k2):
        dis = 0
        kk1 = k1.split("\t")
        kk2 = k2.split("\t")
        for i in range(len(kk1)):
            if (kk1[i] != kk2[i]):
               dis += 1
        return dis
    def __unpickleobj(self, ext):
        f = open(ext, "r")
        t = pickle.load(f)
        f.close()
        return t

    def __pickleobj(self, obj, ext):
        f = open(ext, "w")
        pickle.dump(obj, f)
        f.close()
 
    def loadIdx(self, fname):
        self.keydic = self.__unpickleobj(fname + ".keydic")
        self.seqdic = self.__unpickleobj(fname + ".seqdic")
        self.seqinfo = self.__unpickleobj(fname + ".seqinfo")
        self.seqdata = self.__unpickleobj(fname + ".seqdata")

    def saveIdx(self, fname):
        self.__pickleobj(self.keydic, fname + ".keydic")
        self.__pickleobj(self.seqdic, fname + ".seqdic")
        self.__pickleobj(self.seqinfo, fname + ".seqinfo")
        self.__pickleobj(self.seqdata, fname + ".seqdata")
        

    def addSequence(self, sequence, seqinfo):
        ln = len(sequence)
        seqid = self.seqnumber  
        self.seqnumber += 1
        self.seqinfo.append(seqinfo)
        self.seqdata.append(sequence)
        for i in range(ln - self.windowsize):
            key = "\t".join(sequence[i:(i + self.windowsize)]) 
            self.seqdic.setdefault(key, []).append((seqid, i))
            if (not self.keydic.has_key(key)):
                self.keydic[key] = {}
                self.keydic[key][key] = 0 
                for k in self.keydic.keys(): 
                    if (k != key):
                        dis = self.__keyDistance(k, key)
                        if (dis < 3):
                           self.keydic[key][k] = dis  
                           self.keydic[k][key] = dis  
                           
    def findSequence(self, sequence):
        ln = len(sequence)
        for i in range(ln - self.windowsize):
            key = "\t".join(sequence[i:(i + self.windowsize)]) 
            if (self.keydic.has_key(key)):
                for k in self.keydic[key].keys():
                    dis = self.keydic[key][k]
                    for (id, pos) in self.seqdic[k]:
                        print "pos ", i, self.seqinfo[id], pos, dis
    
def main():
  
    usage = "%prog [options] <datafile.msidx> <datafile.idxed>"
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"

    oparser = optparse.OptionParser(usage=usage, version=version)
    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    idxer = MSIndexer()
    fin = open(args[1], "r")
    while 1:
            line = fin.readline()
            if not line: break;
            line = line.strip("\n\r\t ")
            if (len(line) == 0): continue
            if (line[0]  == ">"):
                cmm = line 
            else:
                seq = line.split(" ")
                print >> sys.stderr, "indxing %s " % cmm
                idxer.addSequence(seq, cmm)
    fin.close()

    # search
    fin = open(args[1], "r")
    for k in range(2):
        line = fin.readline()
        if not line: break;
        line = line.strip("\n\r\t ")
        if (len(line) == 0): continue
        if (line[0]  == ">"):
            cmm = line 
        else:
            seq = line.split(" ")
            idxer.findSequence(seq)
    fin.close()
    idxer.saveIdx(args[2])
    
       
        
    
if __name__ == '__main__':
    main()
