#!/usr/bin/env python

import sys, math, string, optparse, fileinput, struct, pickle

class  MSIndexer:
    def __init__(self):
        self.keydic = {}
        self.seqdic = {}
        self.seqinfo = [] 
        self.seqdata = [] 
        self.seqnumber = 0
        self.keywordlength = 8
        self.keywordcount = 0

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
        self.entry = self.__unpickleobj(fname + ".entry")

    def saveIdx(self, fname):
        self.__pickleobj(self.keydic, fname + ".keydic")
        self.__pickleobj(self.seqdic, fname + ".seqdic")
        self.__pickleobj(self.seqinfo, fname + ".seqinfo")
        self.__pickleobj(self.seqdata, fname + ".seqdata")
        self.__pickleobj(self.entry, fname + ".entry")

    def getkeyfromseq(self, sequence):
        ln = len(sequence)
        keys = []
        for i in range(0, ln - 4 * self.keywordlength, self.keywordlength*2):
            key = "\t".join(sequence[i:(i + 4*self.keywordlength):4]) 
            keys.append((key, i))
            self.keywordcount += 1 
        return keys
        
    def addSequence(self, sequence, seqinfo):
        seqid = self.seqnumber
        self.seqnumber += 1
        self.seqinfo.append(seqinfo)
        self.seqdata.append(sequence)
        keys = self.getkeyfromseq(sequence)
        for k,pos in keys:
            self.seqdic.setdefault(k, []).append((seqid, pos))

    def xcombinations(self, items, n):
        if n==0: yield []
        else:
            for i in xrange(len(items)):
                for cc in self.xcombinations(items[:i]+items[i+1:],n-1):
                    yield [items[i]]+cc

    def similar(self, k1, k2):
        K1 = [int(k) for k in k1.split("\t")]
        K2 = [int(k) for k in k2.split("\t")]
        s = 0.0
        for i in range(self.keywordlength):
            s += max(0, 2 - self.dismatrix[K1[i]][K2[i]])
        if (s > 6):
            return True
        else: 
            return False
            
    def createEntryTable(self, keylength):
        codesize = len(self.dismatrix)
        print >> sys.stderr, "Creating entry table ..., code size %d, keylength %d" %(codesize, keylength)
        codes = ["%d"%c for c in range(codesize)]
        self.entry ={}
        avg = 0.0
        idx = 1
        for c in self.xcombinations(codes, keylength):
            if (int(idx/100) * 100 == idx):
                print >> sys.stderr, "\r %d %.3f" % (idx, avg/idx),
            ekey = "\t".join(c) 
            for k in self.seqdic.keys():
                if self.similar(ekey, k):
                    self.entry.setdefault(ekey, []).append(k)
                    avg += 1
            idx += 1
        print >> sys.stderr, "Done"
        
    def showinfo(self):
        sum = 0.0
        for k in self.seqdic.keys():
            sum += len(self.seqdic[k])
            if (len(self.seqdic[k]) > 2):
                print "\nkey:%s:" % (k)
                for id,pos in self.seqdic[k]:
                    print "%d,%d,%s" %(id, pos, self.seqinfo[id]),
        print "avg length for each k is %.3f" %  (sum / len(self.seqdic.keys()))
        print "%d of uniq keywords, %d of total keywords" % (len(self.seqdic.keys()), self.keywordcount)
         

    def extractKeywords(self):
        if (not self.keydic.has_key(key)):
            self.keydic[key] = {}
            self.keydic[key][key] = 0 
            for k in self.keydic.keys(): 
                if (k != key):
                    dis = self.__keyDistance(k, key)
                    if (dis < 3):
                       self.keydic[key][k] = dis  
                       self.keydic[k][key] = dis  

    def loadcodebook(self, file):
        allcodes = []
        for line in fileinput.input(file):
            line = line.strip("\n\t ")
            fld = [float(t) for t in line.split("\t")]
            allcodes.append(fld)
        self.dismatrix = []
        for x in allcodes:
            dvec = []
            for y in allcodes:
                d = 0.0
                for i in range(len(x)):
                    d += math.fabs(x[i] - y[i]) 
                dvec.append(d)
            self.dismatrix.append(dvec)

    def queryDB(self, sequence):
        #print self.entry.keys()
        cnt = 0
        cntkey = 0
        for en, pos in self.getkeyfromseq(sequence):
            entrylist = self.entry.setdefault(en, []) 
            if (len(entrylist)  > 0):
                print "\nEntry key [%s]: " % (en)
                cntkey += 1
            for k in entrylist:
                print "\nkey [%s] : shapes " %k, 
                for id, pos in self.seqdic[k]:
                    print " (%d, %s) " %( id, self.seqinfo[id] ),
                    cnt += 1
        print  "keycount, shapecount %d %d " %(cntkey, cnt)
    

def main():
  
    usage = "%prog [options] <datafile.msidx> <datafile.idxed>"
    version = "%prog 0.2\nBuild an Indexed MSS database.\nLongbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-k', '--keywordlength', dest = 'keywordlength', type='int',default = 8 , help = 'keyword length, default: 8')
    oparser.add_option('-c', '--codebook', dest = 'codebook', default = None , help = 'codebook file,  default: None')

    (options, args) = oparser.parse_args(sys.argv)
    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    
    idxer = MSIndexer()
    idxer.keywordlength = options.keywordlength
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
    idxer.showinfo()
    if (options.codebook):
        idxer.loadcodebook(options.codebook)
        idxer.createEntryTable(options.keywordlength)
    idxer.saveIdx(args[2])
    
if __name__ == '__main__':
    main()
