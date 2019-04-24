# matplotlib.use("Agg") # Use Agg if you want to draw figs on numerical servers
import matplotlib.pyplot as plt

from numpy import double, arange, floor, zeros, pi, sin, sqrt, tile, transpose, exp, shape
from scipy import fftpack
from mymodule import *
from numpy import max as npmax

class Corr(object):
    """
    Member variables
    ----------------

    Member functions
    ----------------

    Sample
    ------
    """
    def __init__(self, sd, resol, rel=False):
        """
        sd : data of SwarmData class
        resol : determine resolution of raidal direction
        """
        self.sd = sd
        self.pos = sd.pos
        self.u    = sd.u
        self.N = sd.N
        self.timelen = sd.timelen
        
        self.rad = sd.rad
        self.relr = sd.relrad
        self.vel = sd.vel
        self.speed = sd.speed
        
        self.resol = resol
        self.Rmax = max(sd.rad.flatten())*2
        self.relRmax = max(sd.relrad.flatten())*2
        # self.Rmax = npmax(sd.rad, axis=1)*2
        # self.relRmax = npmax(sd.relrad, axis=1)*2
        self.dr = self.Rmax / self.resol
        self.drRel = self.relRmax / self.resol
        self.drs = arange(0, self.Rmax +  self.dr,  self.dr)
        # print self.dr, len(self.drs), self.drs[-1], self.Rmax
        self.drsRel = arange(0, self.relRmax + self.drRel, self.drRel)

        # parameter for smoothed delta function
        sigma = 0.1
        self.s1 = 2.0*sigma*sigma
        self.s2 = 1.0/(sigma*sqrt(2.0*pi))

        if rel == False:
            self.corr()
        else:
            self.corrRel()
    def corr(self):
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)
        self.Cr = zeros([self.timelen,self.resol+1])
        self.Rij = zeros([self.timelen, self.N, self.N])
        self.DotMij = zeros([self.timelen, self.N, self.N])

        for t in range(self.timelen):
            MatX = tile(self.pos[t,Xindxs],(self.N,1))
            MatY = tile(self.pos[t,Yindxs],(self.N,1))
            diffMatX = transpose(MatX) - MatX
            diffMatY = transpose(MatY) - MatY
            self.Rij[t,:,:] = sqrt(diffMatX*diffMatX+diffMatY*diffMatY)

            ux = tile(self.u[t,Xindxs],(self.N,1))
            uy = tile(self.u[t,Yindxs],(self.N,1))
            self.DotMij[t,:,:] = ux*transpose(ux)+uy*transpose(uy)
        
            c0 = sum((self.SmoothDeltaFunc(self.Rij[t,:,:])*self.DotMij[t,:,:]).flatten()) / sum((self.SmoothDeltaFunc(self.Rij[t,:,:])).flatten())
            for i,r in enumerate(self.drs[:(self.resol+1)]):
                self.Cr[t,i] = sum((self.SmoothDeltaFunc(r - self.Rij[t,:,:])*self.DotMij[t,:,:]).flatten()) / (c0*sum((self.SmoothDeltaFunc(r - self.Rij[t,:,:])).flatten()))

    def corrRel(self):
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)
        self.CrRel = zeros([self.timelen,self.resol+1])
        self.RijRel = zeros([self.timelen, self.N, self.N])
        self.DotMijRel = zeros([self.timelen, self.N, self.N])

        for t in range(self.timelen):
            MatX = tile(self.pos[t,Xindxs],(self.N,1))
            MatY = tile(self.pos[t,Yindxs],(self.N,1))
            diffMatX = transpose(MatX) - MatX
            diffMatY = transpose(MatY) - MatY
            self.RijRel[t,:,:] = sqrt(diffMatX*diffMatX+diffMatY*diffMatY)/self.sd.aveR[t]

            ux = tile(self.u[t,Xindxs],(self.N,1))
            uy = tile(self.u[t,Yindxs],(self.N,1))
            self.DotMijRel[t,:,:] = ux*transpose(ux)+uy*transpose(uy)
        
            c0 = sum((self.SmoothDeltaFunc(self.RijRel[t,:,:])*self.DotMijRel[t,:,:]).flatten()) / sum((self.SmoothDeltaFunc(self.RijRel[t,:,:])).flatten())
            for i,r in enumerate(self.drsRel):
                self.CrRel[t,i] = sum((self.SmoothDeltaFunc(r - self.RijRel[t,:,:])*self.DotMijRel[t,:,:]).flatten()) / (c0*sum((self.SmoothDeltaFunc(r - self.RijRel[t,:,:])).flatten()))

    def SmoothDeltaFunc(self, r):
        return exp( - r*r/ self.s1 ) * self.s2
                     

if __name__=='__main__':
    print 'Sample Usage'
