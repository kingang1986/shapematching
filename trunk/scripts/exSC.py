#!/usr/bin/env python

# import the necessary things for OpenCV

import sys, math, optparse, string, fileinput

from numpy import *
from MSS import *


class exSC:

    def __init__(self):
         self.nbins_r = 5
         self.nbins_theta = 12

    def __logspace(self, d1, d2, n):
        sp =  [( 10 **(d1 + k * (d2-d1)/(n - 1)))   for k in range(0, n -1)]
        sp.append(10 ** d2)
        return sp

    def __dist2(self, x, c):
        [ndata, dimx] = x.shape
        [ncentres, dimc] = c.shape
        result = zeros((ndata, ncentres))
        for i in range(ndata):
            for j in range(ncentres):
                ddata = x[i, :] - c[j, :]
                result[i,j] = sqrt(sum(ddata * ddata))
        return result

    def ExtractPartSC(self, Bsamp, Psamp):
        [nsamp, ndim] = Bsamp.shape
        [npart, ndim] = Psamp.shape
        nbins_theta = self.nbins_theta
        nbins_r = self.nbins_r
        r_inner = 1.0/8
        r_outer = 2.0
        r_array = self.__dist2(Bsamp, Psamp)
        #dy = Bsamp[:,1]*ones((1,nsamp)) - ones((nsamp,1)) * (Bsamp[:,1].transpose())
        dy = dot(ones((nsamp,1)) , (Psamp[:,1:2].transpose())) - dot(Bsamp[:,1:2], ones((1, npart))) 
        dx = dot(ones((nsamp,1)) , (Psamp[:,0:1].transpose())) - dot(Bsamp[:,0:1], ones((1, npart)))
        #theta_array = arctan2(Bsamp[:,1]*ones((1,nsamp)) - ones((nsamp,1)) * Bsamp[:,0], Bsamp[:, 0].transpose() * ones((1,nsamp)) - ones((nsamp,1)) * Bsamp[:, 0]).transpose() 
        theta_array = arctan2(dy, dx)
        mean_dist = r_array.mean()
        r_array_n = r_array / mean_dist
        r_bin_edges = self.__logspace(math.log10(r_inner), math.log10(r_outer), 5)
        r_array_q = zeros((nsamp,npart), dtype=int)
        for m in range(self.nbins_r):
           r_array_q = r_array_q + (r_array_n < r_bin_edges[m])
        fz = r_array_q > 0  #  flag all points inside outer boundary
        
        # put all angles in [0,2pi) range
        theta_array_2 = theta_array + 2* math.pi * (theta_array < 0) 
        
        # quantize to a fixed set of angles (bin edges lie on 0,(2*pi)/k,...2*pi
        theta_array_q = 1 + floor(theta_array_2 /(2 * math.pi / nbins_theta))
        nbins = nbins_theta*nbins_r
        BH = zeros((nsamp,nbins))
        for n in range(nsamp):
            sn = zeros((nbins_r, nbins_theta))
            for k in range(npart):
                if (fz[n, k]):
                    sn[r_array_q[n, k] - 1, theta_array_q[n, k] - 1] += 1
            BH[n,:] = sn.reshape(nbins)
        return BH
        
    def ExtractSC(self, Bsamp):
        [nsamp, ndim] = Bsamp.shape
        nbins_theta = self.nbins_theta
        nbins_r = self.nbins_r
        r_inner = 1.0/8
        r_outer = 2.0
        r_array = self.__dist2(Bsamp, Bsamp)
        #dy = Bsamp[:,1]*ones((1,nsamp)) - ones((nsamp,1)) * (Bsamp[:,1].transpose())
        dy = dot(ones((nsamp,1)) , (Bsamp[:,1:2].transpose())) - dot(Bsamp[:,1:2], ones((1, nsamp))) 
        dx = dot(ones((nsamp,1)) , (Bsamp[:,0:1].transpose())) - dot(Bsamp[:,0:1], ones((1, nsamp)))
        #theta_array = arctan2(Bsamp[:,1]*ones((1,nsamp)) - ones((nsamp,1)) * Bsamp[:,0], Bsamp[:, 0].transpose() * ones((1,nsamp)) - ones((nsamp,1)) * Bsamp[:, 0]).transpose() 
        theta_array = arctan2(dy, dx)
        mean_dist = r_array.mean()
        r_array_n = r_array / mean_dist
        r_bin_edges = self.__logspace(math.log10(r_inner), math.log10(r_outer), 5)
        r_array_q = zeros((nsamp,nsamp), dtype=int)
        for m in range(self.nbins_r):
           r_array_q = r_array_q + (r_array_n < r_bin_edges[m])
        fz = r_array_q > 0  #  flag all points inside outer boundary
        
        # put all angles in [0,2pi) range
        theta_array_2 = theta_array + 2* math.pi * (theta_array < 0) 
        
        # quantize to a fixed set of angles (bin edges lie on 0,(2*pi)/k,...2*pi
        theta_array_q = 1 + floor(theta_array_2 /(2 * math.pi / nbins_theta))
        nbins = nbins_theta*nbins_r
        BH = zeros((nsamp,nbins))
        for n in range(nsamp):
            sn = zeros((nbins_r, nbins_theta))
            for k in range(nsamp):
                if (fz[n, k]):
                    sn[r_array_q[n, k] - 1, theta_array_q[n, k] - 1] += 1
            BH[n,:] = sn.reshape(nbins)
        return BH


    def ExtractFeature(self, keypoints, targetlist):
        nsamp = len(keypoints)
        Bsamp = ones((nsamp, 2))
        tsamp = nsamp
        for k in range(len(keypoints)):
            Bsamp[k, 0] = keypoints[k].x
            Bsamp[k, 1] = keypoints[k].y
        if (targetlist):
            Psamp = Bsamp[targetlist]
            BH = self.ExtractPartSC(Bsamp, Psamp)
            tsamp = len(targetlist) 
        else:
            BH = self.ExtractSC(Bsamp)

        varname = ["fsc%2d"%(s) for s in range(self.nbins_theta * self.nbins_r)]
        for s in range(self.nbins_theta * self.nbins_r):
            for k in range(len(keypoints)):
                keypoints[k].addvalue(varname[s], BH[k, s]/tsamp)
        varname = ["gsc%2d"%(s) for s in range(self.nbins_theta * self.nbins_r)]
        for s in range(self.nbins_theta * self.nbins_r):
            for k in range(len(keypoints)):
                keypoints[k].addvalue(varname[s], BH[k, s]/tsamp)
 
def main():
    usage = "%s [options]  <MSS File> <MSSout>" % (sys.argv[0])
    version = "%prog 0.2\n Longbin Chen, longbinc@yahoo.com"
    oparser = optparse.OptionParser(usage=usage, version=version)
    oparser.add_option('-t', '--head', dest = 'head', type='int',default = None, help = 'use first h as target')

    (options, args) = oparser.parse_args(sys.argv)

    if len(args) != 3:
        oparser.parse_args([sys.argv[0], "--help"])
        sys.exit(1)
    sc = exSC()
    mss = MSS()
    mss.load(args[1])
    keypoints = mss.getAllPoints()
    if (options.head == None):
        hl = len(keypoints)
    else:
        hl = options.head
    sc.ExtractFeature(keypoints, range(hl))
    mss.save(args[2])

if __name__ == '__main__':
    # load the image gived on the command line
    main()

