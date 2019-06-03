from struct import unpack
from sys import exit

from numpy import array, double, arange, float64, zeros, sqrt, loadtxt, tile, transpose, fromfile, mean, sin, cos, arctan2, sum, add, shape
from mymodule.LyapunovAnalysis import *
import mymodule.linearft as lf

"""
class StaticData(object):   class for -SwarmState file
class LyapunovData(object): class for -LyapunovSpectra file
class IndLyapData(object):  class for -IndLyapunovSpectra file
class SwarmData(object):    class for -State file (binary file)
class ForceTS(object):    class for -State file (binary file)
"""

class StaticData(object):
    """
    class for SwarmState file
    class StaticData(object):   class for -SwarmState file
    """
    def __init__(self, filename):
        """
        Parameter
        ---------
        filename : filename of SwarmState file
                   if filename is empyt, this class is assumed for statistical amount
        """
        self.filename = filename
        if not self.filename:
            self.N            = 0
            self.Lmax         = 0.0
            self.Vcom         = 0.0
            self.Ang          = 0.0
            self.AbsAng       = 0.0
            self.NormedAng    = 0.0
            self.AbsNormedAng = 0.0
            self.Polarity     = 0.0
            self.phi          = 0.0
            self.Ep           = 0.0
            self.Ek           = 0.0
            self.Size         = 0.0
        else:
            self.N = int(filename[(filename.find('Num')+3):filename.find('v0')])
            self.errorcheck(self.filename)
        
            self.Lmax         = self.data[8]*100.0
            self.Vcom         = self.data[10]*100.0
            self.Ang          = self.data[12]*100.0
            self.AbsAng       = self.data[30]*100.0
            self.NormedAng    = abs(self.data[32]*100.0)
            self.AbsNormedAng = self.data[34]*100.0
            self.Polarity     = self.data[24]*100.0
            self.phi          = self.data[14]*100.0
            self.Ep           = self.data[26]*100.0
            self.Ek           = self.data[28]*100.0
            self.Size         = self.data[20]*100.0

    def errorcheck(self, filename):
        if filename.find('SwarmState')==-1:
            print("This file is not SwarmState file ", filename)
            exit(1)

        try:
            self.data=loadtxt(filename,skiprows=1)
        except IOError:
            print("Fail to load file ", filename)
            exit(1)
    
    def __add__(self, rhs):
        sd = StaticData('')
        sd.Lmax         = self.Lmax         + rhs.Lmax
        sd.Vcom         = self.Vcom         + rhs.Vcom
        sd.Ang          = self.Ang          + rhs.Ang
        sd.AbsAng       = self.AbsAng       + rhs.AbsAng
        sd.NormedAng    = self.NormedAng    + rhs.NormedAng
        sd.AbsNormedAng = self.AbsNormedAng + rhs.AbsNormedAng
        sd.Polarity     = self.Polarity     + rhs.Polarity
        sd.phi          = self.phi          + rhs.phi
        sd.Ep           = self.Ep           + rhs.Ep
        sd.Ek           = self.Ek           + rhs.Ek
        sd.Size         = self.Size         + rhs.Size

        return sd

    def __sub__(self, rhs):
        sd = StaticData('')
        sd.Lmax         = self.Lmax         - rhs.Lmax
        sd.Vcom         = self.Vcom         - rhs.Vcom
        sd.Ang          = self.Ang          - rhs.Ang
        sd.AbsAng       = self.AbsAng       - rhs.AbsAng
        sd.NormedAng    = self.NormedAng    - rhs.NormedAng
        sd.AbsNormedAng = self.AbsNormedAng - rhs.AbsNormedAng
        sd.Polarity     = self.Polarity     - rhs.Polarity
        sd.phi          = self.phi          - rhs.phi
        sd.Ep           = self.Ep           - rhs.Ep
        sd.Ek           = self.Ek           - rhs.Ek
        sd.Size         = self.Size         - rhs.Size

        return sd

    def __div__(self, rhs):
        sd = StaticData('')
        sd.Lmax         = self.Lmax         / rhs
        sd.Vcom         = self.Vcom         / rhs
        sd.Ang          = self.Ang          / rhs
        sd.AbsAng       = self.AbsAng       / rhs
        sd.NormedAng    = self.NormedAng    / rhs
        sd.AbsNormedAng = self.AbsNormedAng / rhs
        sd.Polarity     = self.Polarity     / rhs
        sd.phi          = self.phi          / rhs
        sd.Ep           = self.Ep           / rhs
        sd.Ek           = self.Ek           / rhs
        sd.Size         = self.Size         / rhs

        return sd

    def __mul__(self, rhs):
        sd = StaticData('')
        sd.Lmax         = self.Lmax         * rhs
        sd.Vcom         = self.Vcom         * rhs
        sd.Ang          = self.Ang          * rhs
        sd.AbsAng       = self.AbsAng       * rhs
        sd.NormedAng    = self.NormedAng    * rhs
        sd.AbsNormedAng = self.AbsNormedAng * rhs
        sd.Polarity     = self.Polarity     * rhs
        sd.phi          = self.phi          * rhs
        sd.Ep           = self.Ep           * rhs
        sd.Ek           = self.Ek           * rhs
        sd.Size         = self.Size         * rhs

        return sd

    def __iadd__(self, rhs):
        self.Lmax         += rhs.Lmax
        self.Vcom         += rhs.Vcom
        self.Ang          += rhs.Ang
        self.AbsAng       += rhs.AbsAng
        self.NormedAng    += rhs.NormedAng
        self.AbsNormedAng += rhs.AbsNormedAng
        self.Polarity     += rhs.Polarity
        self.phi          += rhs.phi
        self.Ep           += rhs.Ep
        self.Ek           += rhs.Ek
        self.Size         += rhs.Size

        return self

    def __isub__(self, rhs):
        self.Lmax         -= rhs.Lmax
        self.Vcom         -= rhs.Vcom
        self.Ang          -= rhs.Ang
        self.AbsAng       -= rhs.AbsAng
        self.NormedAng    -= rhs.NormedAng
        self.AbsNormedAng -= rhs.AbsNormedAng
        self.Polarity     -= rhs.Polarity
        self.phi          -= rhs.phi
        self.Ep           -= rhs.Ep
        self.Ek           -= rhs.Ek
        self.Size         -= rhs.Size

        return self

    def __idiv__(self, rhs):
        self.Lmax         /= rhs
        self.Vcom         /= rhs
        self.Ang          /= rhs
        self.AbsAng       /= rhs
        self.NormedAng    /= rhs
        self.AbsNormedAng /= rhs
        self.Polarity     /= rhs
        self.phi          /= rhs
        self.Ep           /= rhs
        self.Ek           /= rhs
        self.Size         /= rhs

        return self
            
    def __imul__(self, rhs):
        self.Lmax         *= rhs
        self.Vcom         *= rhs
        self.Ang          *= rhs
        self.AbsAng       *= rhs
        self.NormedAng    *= rhs
        self.AbsNormedAng *= rhs
        self.Polarity     *= rhs
        self.phi          *= rhs
        self.Ep           *= rhs
        self.Ek           *= rhs
        self.Size         *= rhs

        return self

class LyapunovData(object):
    """
    class for LyapunovSpectra file
    class LyapunovData(object): class for -LyapunovSpectra file
    """
    def __init__(self, filename, dim=4):
        self.filename = filename
        if not self.filename:
            self.LES = array([])
            self.DL = 0.0
            self.pLEs = 0.0
            self.nLEs = 0.0
            self.pos_LEs = 0.0
            self.neg_LEs = 0.0
            self.np_LEs = 0.0
        else:
            self.N = int(filename[(filename.find('Num')+3):filename.find('v0')])
            self.errorcheck(self.filename)

            tmp = len(self.data[:,0])-dim*self.N
            self.LES = array(self.data[tmp:,1])
            self.DL = CalcDky(self.LES)/double(dim*self.N)
            self.pLEs = sum(self.LES[self.LES>0])
            self.nLEs = abs(sum(self.LES[self.LES<0]))
            self.pos_LEs = self.pLEs/double(self.N)
            self.neg_LEs = self.nLEs/double(self.N)
            self.np_LEs = self.pLEs/self.nLEs
        
    def errorcheck(self, filename):
        if filename.find('LyapunovSpectra')==-1:
            print("This file is not SwarmState file ", filename)
            exit(1)

        try:
            self.data=loadtxt(filename)
        except IOError:
            print("Fail to load file ", filename)
            exit(1)

    def __add__(self, rhs):
        ld = LyapunovData('')
        ld.LES      = self.LES     + rhs.LES
        ld.DL       = self.DL      + rhs.DL
        ld.pLEs     = self.pLEs    + rhs.pLEs
        ld.nLEs     = self.nLEs    + rhs.nLEs
        ld.pos_LEs  = self.pos_LEs + rhs.pos_LEs
        ld.neg_LEs  = self.neg_LEs + rhs.neg_LEs
        ld.np_LEs   = self.np_LEs  + rhs.np_LEs
        return ld

    def __sub__(self, rhs):
        ld = LyapunovData('')
        ld.LES      = self.LES     - rhs.LES
        ld.DL       = self.DL      - rhs.DL
        ld.pLEs     = self.pLEs    - rhs.pLEs
        ld.nLEs     = self.nLEs    - rhs.nLEs
        ld.pos_LEs  = self.pos_LEs - rhs.pos_LEs
        ld.neg_LEs  = self.neg_LEs - rhs.neg_LEs
        ld.np_LEs   = self.np_LEs  - rhs.np_LEs
        return ld

    def __div__(self, rhs):
        ld = LyapunovData('')
        ld.LES      = self.LES     / rhs
        ld.DL       = self.DL      / rhs
        ld.pLEs     = self.pLEs    / rhs
        ld.nLEs     = self.nLEs    / rhs
        ld.pos_LEs  = self.pos_LEs / rhs
        ld.neg_LEs  = self.neg_LEs / rhs
        ld.np_LEs   = self.np_LEs  / rhs
        return ld

    def __mul__(self, rhs):
        ld = LyapunovData('')
        ld.LES      = self.LES     * rhs
        ld.DL       = self.DL      * rhs
        ld.pLEs     = self.pLEs    * rhs
        ld.nLEs     = self.nLEs    * rhs
        ld.pos_LEs  = self.pos_LEs * rhs
        ld.neg_LEs  = self.neg_LEs * rhs
        ld.np_LEs   = self.np_LEs  * rhs
        return ld

    def __iadd__(self, rhs):
        self.LES     += rhs.LES
        self.DL      += rhs.DL
        self.pLEs    += rhs.pLEs
        self.nLEs    += rhs.nLEs
        self.pos_LEs += rhs.pos_LEs
        self.neg_LEs += rhs.neg_LEs
        self.np_LEs  += rhs.np_LEs
        return self

    def __isub__(self, rhs):
        self.LES     -= rhs.LES
        self.DL      -= rhs.DL
        self.pLEs    -= rhs.pLEs
        self.nLEs    -= rhs.nLEs
        self.pos_LEs -= rhs.pos_LEs
        self.neg_LEs -= rhs.neg_LEs
        self.np_LEs  -= rhs.np_LEs
        return self

    def __idiv__(self, rhs):
        self.LES     /= rhs
        self.DL      /= rhs
        self.pLEs    /= rhs
        self.nLEs    /= rhs
        self.pos_LEs /= rhs
        self.neg_LEs /= rhs
        self.np_LEs  /= rhs
        return self

    def __imul__(self, rhs):
        self.LES     *= rhs
        self.DL      *= rhs
        self.pLEs    *= rhs
        self.nLEs    *= rhs
        self.pos_LEs *= rhs
        self.neg_LEs *= rhs
        self.np_LEs  *= rhs
        return self

class IndLyapData(object):
    """
    class for IndLyapunovSpectra file
    class IndLyapData(object):  class for -IndLyapunovSpectra file
    """
    def __init__(self, filename, dim=4):
        self.filename = filename
        if not self.filename:
            self.LES = array([])
            self.DL = array([])
            self.pLEs = array([])
            self.nLEs = array([])
            self.np_LEs = array([])
        else:
            self.N = int(filename[(filename.find('Num')+3):filename.find('v0')])
            self.errorcheck(self.filename)

            self.LES = array(self.data)
            self.DL = array([ CalcDky(self.LES[i,:])/double(dim) for i in range(self.N) ])
            self.pLEs = array([ sum(self.LES[i][self.LES[i,:]>0]) for i in range(self.N) ])
            self.nLEs = array([ abs(sum(self.LES[i][self.LES[i,:]<0])) for i in range(self.N) ])
            self.np_LEs = array([ self.pLEs[i]/self.nLEs[i] for i in range(self.N) ])
            
    def errorcheck(self, filename):
        if filename.find('IndLyapunovSpectra')==-1:
            print("This file is not IndLyapunovSpectra file ", filename)
            exit(1)

        try:
            self.data=loadtxt(filename, skiprows=1) # first row is just index [0, 1, 2, 3]
        except IOError:
            print("Fail to load file ", filename)
            exit(1)
    def __add__(self, rhs):
        ld = IndLyapData('')
        ld.LES      = self.LES     + rhs.LES
        ld.DL       = self.DL      + rhs.DL
        ld.pLEs     = self.pLEs    + rhs.pLEs
        ld.nLEs     = self.nLEs    + rhs.nLEs
        ld.np_LEs   = self.np_LEs  + rhs.np_LEs
        return ld

    def __sub__(self, rhs):
        ld = IndLyapData('')
        ld.LES      = self.LES     - rhs.LES
        ld.DL       = self.DL      - rhs.DL
        ld.pLEs     = self.pLEs    - rhs.pLEs
        ld.nLEs     = self.nLEs    - rhs.nLEs
        ld.np_LEs   = self.np_LEs  - rhs.np_LEs
        return ld

    def __div__(self, rhs):
        ld = IndLyapData('')
        ld.LES      = self.LES     / rhs
        ld.DL       = self.DL      / rhs
        ld.pLEs     = self.pLEs    / rhs
        ld.nLEs     = self.nLEs    / rhs
        ld.np_LEs   = self.np_LEs  / rhs
        return ld

    def __mul__(self, rhs):
        ld = IndLyapData('')
        ld.LES      = self.LES     * rhs
        ld.DL       = self.DL      * rhs
        ld.pLEs     = self.pLEs    * rhs
        ld.nLEs     = self.nLEs    * rhs
        ld.np_LEs   = self.np_LEs  * rhs
        return ld

    def __iadd__(self, rhs):
        self.LES     += rhs.LES
        self.DL      += rhs.DL
        self.pLEs    += rhs.pLEs
        self.nLEs    += rhs.nLEs
        self.np_LEs  += rhs.np_LEs
        return self

    def __isub__(self, rhs):
        self.LES     -= rhs.LES
        self.DL      -= rhs.DL
        self.pLEs    -= rhs.pLEs
        self.nLEs    -= rhs.nLEs
        self.np_LEs  -= rhs.np_LEs
        return self

    def __idiv__(self, rhs):
        self.LES     /= rhs
        self.DL      /= rhs
        self.pLEs    /= rhs
        self.nLEs    /= rhs
        self.np_LEs  /= rhs
        return self

    def __imul__(self, rhs):
        self.LES     *= rhs
        self.DL      *= rhs
        self.pLEs    *= rhs
        self.nLEs    *= rhs
        self.np_LEs  *= rhs
        return self


class SwarmData(object):
    """
    Class for binary file of the State data
    class SwarmData(object):    class for -State file (binary file)
    """
    def __init__(self, filename):
        self.filename = filename
        self.N = int(filename[(filename.find('Num')+3):filename.find('v0')])

        self.LoadBinaryState()
        
        self.CalcSwarmStatics()
        self.CalcRelativeCoordinate()

    def LoadBinaryState(self):
        """
        Prameters
        ---------
        f : file object you want to load. f is NOT FILENAME
        N : Number of particles
        
        Returns
        -------
        timelen : duration of saved data
        pos     : global position data
        vel     : global velocity data
        """
        PhaseDim=4*self.N
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)

        with open(self.filename, 'rb') as f:
            sys_info = unpack('<ii',f.read(8))
            tmp = fromfile(f,dtype='<d')
            self.timelen = len(tmp)/(4*sys_info[0])
            print(len(tmp),sys_info)
            tmp = tmp.reshape(self.timelen, 4*sys_info[0])
            self.ps_data = tmp
        
            # self.timelen = len(tmp[:,0])
        
            self.pos = zeros([self.timelen,2*self.N])
            self.vel = zeros([self.timelen,2*self.N])
        
            self.pos[:,Xindxs] = tmp[:,0:PhaseDim:4]
            self.pos[:,Yindxs] = tmp[:,1:PhaseDim:4]
        
            self.vel[:,Xindxs] = tmp[:,2:PhaseDim:4]
            self.vel[:,Yindxs] = tmp[:,3:PhaseDim:4]
            # self.virial = sum( self.pos[:,Xindxs]*self.vel[:,Xindxs] + self.pos[:,Yindxs]*self.vel[:,Yindxs], axis=1 )

    def CalcSwarmStatics(self):
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)
        
        self.Rcom = zeros([self.timelen,2])
        self.Vcom = zeros([self.timelen,2])
        
        self.Rcom[:,0] = mean(self.pos[:,Xindxs], axis=1, dtype=float64)
        self.Rcom[:,1] = mean(self.pos[:,Yindxs], axis=1, dtype=float64)
        self.Vcom[:,0] = mean(self.vel[:,Xindxs], axis=1, dtype=float64)
        self.Vcom[:,1] = mean(self.vel[:,Yindxs], axis=1, dtype=float64)

        self.Vabs = sqrt(self.Vcom[:,0]**2 + self.Vcom[:,1]**2)
        self.Rsq  = self.Rcom[:,0]**2 + self.Rcom[:,1]**2
        self.Thetacom = arctan2(self.Vcom[:,0],self.Vcom[:,1])
        
    def CalcRelativeCoordinate(self):
        """
        Parameters
        ----------
        timelen : duration of saved data
        N       : Number of particles
        pos     : global position
        vel     : global velocity
        
        Returns
        -------
        r : Relative position ( substract the position vector of center of mass )
        u : Relative velocity ( substract the velocity of center of mass )
        w : Coordinated velocity ( substract the velocity of center of mass and rotate )
        """

        ## Index dammies
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)

        ## Initialization
        self.r = zeros([self.timelen,2*self.N])
        self.u = zeros([self.timelen,2*self.N])
        self.brcd = zeros([self.timelen,2*self.N])
        self.w = zeros([self.timelen,2*self.N])
        
        self.rad   = zeros([self.timelen,self.N])
        self.theta = zeros([self.timelen,self.N])

        ## cos matrix and sin matrix
        CosMat = transpose(tile(cos(self.Thetacom),(self.N,1)))
        SinMat = transpose(tile(sin(self.Thetacom),(self.N,1)))

        ## Relative position ( substract the position vector of center of mass )
        self.r[:,Xindxs] = self.pos[:,Xindxs] - transpose(tile(self.Rcom[:,0],[self.N,1]))
        self.r[:,Yindxs] = self.pos[:,Yindxs] - transpose(tile(self.Rcom[:,1],[self.N,1]))

        ## Relative velocity ( substract the velocity of center of mass )
        self.u[:,Xindxs] = self.vel[:,Xindxs] - transpose(tile(self.Vcom[:,0],[self.N,1]))
        self.u[:,Yindxs] = self.vel[:,Yindxs] - transpose(tile(self.Vcom[:,1],[self.N,1]))

        ## Coordinated position ( substract the position vector of center of mass and rotate )
        self.brcd[:,Xindxs] =  CosMat * self.r[:,Xindxs] + SinMat * self.r[:,Yindxs]
        self.brcd[:,Yindxs] = -SinMat * self.r[:,Xindxs] + CosMat * self.r[:,Yindxs]

        ## Coordinated velocity ( substract the velocity of center of mass and rotate )
        self.w[:,Xindxs] =  CosMat * self.u[:,Xindxs] + SinMat * self.u[:,Yindxs]
        self.w[:,Yindxs] = -SinMat * self.u[:,Xindxs] + CosMat * self.u[:,Yindxs]

        ## Angular Momentum against to the center of mass
        self.angmom = self.r[:,Xindxs]*self.vel[:,Yindxs] - self.r[:,Yindxs]*self.vel[:,Xindxs]

        ## Individual speed
        self.speed  = sqrt(self.vel[:,Xindxs]*self.vel[:,Xindxs] + self.vel[:,Yindxs]*self.vel[:,Yindxs])

        ## individual relative velocity
        self.relsp  = sqrt(self.u[:,Xindxs]*self.u[:,Xindxs] + self.u[:,Yindxs]*self.u[:,Yindxs])

        ## Individual radius againt to the center of mass
        self.rad   = sqrt(self.r[:,Xindxs] * self.r[:,Xindxs] + self.r[:,Yindxs] * self.r[:,Yindxs])
        self.aveR = mean(self.rad,axis=1)
        meanrad = tile(self.aveR,(self.N,1)).transpose()
        self.relrad= self.rad/meanrad

        ## Individual angle to the X axis
        self.theta = arctan2(self.r[:,Xindxs],self.r[:,Yindxs])

        ## virial of the system
        self.virial = sum( self.r[:,Xindxs]*self.vel[:,Xindxs] + self.r[:,Yindxs]*self.vel[:,Yindxs], axis=1 )

        ## normalized 
        self.Lnormed = self.angmom/(self.speed*sqrt(self.r[:,Xindxs]**2+self.r[:,Yindxs]**2))
        
    def Average_Coeff_Relative(self):
        ps = [ self.Coeff_RelativeDiffRot(i,j) for i in range(self.N) for j in range(self.N) if i != j ]
        
        return mean(ps,axis=0)

    def Coeff_RelativeDiffRot(self,i,j):
        if i == j:
            print("Error : i and j must be different")
            exit(0)

        init_dist = sqrt( (self.r[0,2*i]-self.r[0,2*j])**2 + (self.r[0,2*i+1]-self.r[0,2*j+1])**2 )
        ref_dist  = add.accumulate( tile(init_dist,self.timelen) )
        dist = sqrt( (self.r[:,2*i]-self.r[:,2*j])**2 + (self.r[:,2*i+1]-self.r[:,2*j+1])**2 )
        rot  = arctan2(self.r[:,2*i]-self.r[:,2*j], self.r[:,2*i+1]-self.r[:,2*j+1])
        diff_ij       = add.accumulate(dist) - ref_dist
        diff_ij_abs   = abs( add.accumulate(dist) - ref_dist )
        diff_ij_ratio = add.accumulate(dist)/ref_dist
        rot_ij        = add.accumulate(rot)
        #################################################################
        # self.dist_rel4init = dist/init_dist
        
        # ftdata=lf.LeastSquareOptimize(xdata, self.dist_rel4init, 'linear')
        # ftdata.optimize(1.0,0.0)
        # p_d = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        #################################################################

        xdata= arange(0,len(dist))

        ftdata=lf.LeastSquareOptimize(xdata, diff_ij, 'linear')
        ftdata.optimize(1.0,0.0)
        p = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        ftdata=lf.LeastSquareOptimize(xdata, diff_ij_abs, 'linear')
        ftdata.optimize(1.0,0.0)
        p_abs = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        ftdata=lf.LeastSquareOptimize(xdata, diff_ij_ratio, 'linear')
        ftdata.optimize(1.0,0.0)
        p_ratio = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        ftdata=lf.LeastSquareOptimize(xdata, rot_ij, 'linear')
        ftdata.optimize(1.0,0.0)
        p_rot = ftdata.GetOptParams()

        # self.rel_diff_coeff = [p[0], p_abs[0], p_ratio[0], p_rot[0]]
        return p[0], p_abs[0], p_ratio[0], p_rot[0]

        # ts = fac.TimeSeries(dist,1)
        # ts.ShowPowerSpectorum()

        # plt.close()
        # [hist,bins] = histogram(dist,bins=100)
        # plt.bar(bins[0:len(hist)],hist,width=bins[1]-bins[0])
        # plt.show()

class ForceTS(object):
    """
    Class for binary file of the State data
    class ForceTS(object):    class for -State file (binary file)
    """
    def __init__(self, filename):
        self.filename = filename
        self.N = int(filename[(filename.find('Num')+3):filename.find('v0')])
        self.ifs = open(filename, 'rb')
        self.LoadBinaryState()

        self.f.close()
        
    def LoadBinaryState(self):
        """
        Prameters
        ---------
        ifs : file object you want to load. f is NOT FILENAME
        N : Number of particles
        
        Returns
        -------
        timelen : duration of saved data
        f      : global position data
        df     : global velocity data
        """
        items = 5
        PhaseDim = items*self.N
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)
        
        tmp = fromfile(self.ifs,dtype='<d')
        self.timelen = shape(tmp)[0]/(items*self.N)
        tmp = tmp.reshape(self.timelen,items*self.N);

        self.f = zeros([self.timelen,2*self.N])
        self.df = zeros([self.timelen,2*self.N])
        self.kappa = zeros([self.timelen,self.N])
        
        self.f[:,Xindxs] = tmp[:,0:PhaseDim:items]
        self.f[:,Yindxs] = tmp[:,1:PhaseDim:items]
        self.fabs = sqrt( self.f[:,Xindxs]**2 + self.f[:,Yindxs]**2)
        
        self.df[:,Xindxs] = tmp[:,2:PhaseDim:items]
        self.df[:,Yindxs] = tmp[:,3:PhaseDim:items]
        self.dfabs = sqrt( self.df[:,Xindxs]**2 + self.df[:,Yindxs]**2)

        self.kappa[:,arange(self.N)] = tmp[:,4:PhaseDim:items]
        # self.virial = sum( self.pos[:,Xindxs]*self.vel[:,Xindxs] + self.pos[:,Yindxs]*self.vel[:,Yindxs], axis=1 )

class FisheryData(object):
    """
    Class for binary file of the State data
    class SwarmData(object):    class for -State file (binary file)
    """
    def __init__(self, filename):
        self.filename = filename
        self.N = int(filename[(filename.find('N')+1):filename.find('T')])
        self.f = open(filename, 'rb')
        self.namecore = filename[:filename.find('.data')]

        self.LoadBinaryState()
        
        self.CalcSwarmStatics()
        self.CalcRelativeCoordinate()

        self.f.close()

    def LoadBinaryState(self):
        """
        Prameters
        ---------
        f : file object you want to load. f is NOT FILENAME
        N : Number of particles
        
        Returns
        -------
        timelen : duration of saved data
        pos     : global position data
        vel     : global velocity data
        """
        PhaseDim=4*self.N
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)
        
        sys_info = unpack('<ii',self.f.read(8))
        
        tmp = fromfile(self.f,dtype='<d')
        tmp = tmp.reshape(sys_info[1]+1,4*sys_info[0]);
        self.ps_data = tmp
        
        self.timelen = len(tmp[:,0])
        
        self.pos = zeros([self.timelen,2*self.N])
        self.vel = zeros([self.timelen,2*self.N])
        
        self.pos[:,Xindxs] = tmp[:,0:PhaseDim:4]
        self.pos[:,Yindxs] = tmp[:,1:PhaseDim:4]
        
        self.vel[:,Xindxs] = tmp[:,2:PhaseDim:4]
        self.vel[:,Yindxs] = tmp[:,3:PhaseDim:4]
        # self.virial = sum( self.pos[:,Xindxs]*self.vel[:,Xindxs] + self.pos[:,Yindxs]*self.vel[:,Yindxs], axis=1 )

    def CalcSwarmStatics(self):
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)
        
        self.Rcom = zeros([self.timelen,2])
        self.Vcom = zeros([self.timelen,2])
        
        self.Rcom[:,0] = mean(self.pos[:,Xindxs], axis=1, dtype=float64)
        self.Rcom[:,1] = mean(self.pos[:,Yindxs], axis=1, dtype=float64)
        self.Vcom[:,0] = mean(self.vel[:,Xindxs], axis=1, dtype=float64)
        self.Vcom[:,1] = mean(self.vel[:,Yindxs], axis=1, dtype=float64)

        self.Vabs = sqrt(self.Vcom[:,0]**2 + self.Vcom[:,1]**2)
        self.Rsq  = self.Rcom[:,0]**2 + self.Rcom[:,1]**2
        self.Thetacom = arctan2(self.Vcom[:,0],self.Vcom[:,1])
        
    def CalcRelativeCoordinate(self):
        """
        Parameters
        ----------
        timelen : duration of saved data
        N       : Number of particles
        pos     : global position
        vel     : global velocity
        
        Returns
        -------
        r : Relative position ( substract the position vector of center of mass )
        u : Relative velocity ( substract the velocity of center of mass )
        w : Coordinated velocity ( substract the velocity of center of mass and rotate )
        """

        ## Index dammies
        Xindxs = arange(0,2*self.N,2)
        Yindxs = arange(1,2*self.N,2)

        ## Initialization
        self.r = zeros([self.timelen,2*self.N])
        self.u = zeros([self.timelen,2*self.N])
        self.brcd = zeros([self.timelen,2*self.N])
        self.w = zeros([self.timelen,2*self.N])
        
        self.rad   = zeros([self.timelen,self.N])
        self.theta = zeros([self.timelen,self.N])

        ## cos matrix and sin matrix
        CosMat = transpose(tile(cos(self.Thetacom),(self.N,1)))
        SinMat = transpose(tile(sin(self.Thetacom),(self.N,1)))

        ## Relative position ( substract the position vector of center of mass )
        self.r[:,Xindxs] = self.pos[:,Xindxs] - transpose(tile(self.Rcom[:,0],[self.N,1]))
        self.r[:,Yindxs] = self.pos[:,Yindxs] - transpose(tile(self.Rcom[:,1],[self.N,1]))

        ## Relative velocity ( substract the velocity of center of mass )
        self.u[:,Xindxs] = self.vel[:,Xindxs] - transpose(tile(self.Vcom[:,0],[self.N,1]))
        self.u[:,Yindxs] = self.vel[:,Yindxs] - transpose(tile(self.Vcom[:,1],[self.N,1]))

        ## Coordinated position ( substract the position vector of center of mass and rotate )
        self.brcd[:,Xindxs] =  CosMat * self.r[:,Xindxs] + SinMat * self.r[:,Yindxs]
        self.brcd[:,Yindxs] = -SinMat * self.r[:,Xindxs] + CosMat * self.r[:,Yindxs]

        ## Coordinated velocity ( substract the velocity of center of mass and rotate )
        self.w[:,Xindxs] =  CosMat * self.u[:,Xindxs] + SinMat * self.u[:,Yindxs]
        self.w[:,Yindxs] = -SinMat * self.u[:,Xindxs] + CosMat * self.u[:,Yindxs]

        ## Angular Momentum against to the center of mass
        self.angmom = self.r[:,Xindxs]*self.vel[:,Yindxs] - self.r[:,Yindxs]*self.vel[:,Xindxs]

        ## Individual speed
        self.speed  = sqrt(self.vel[:,Xindxs]*self.vel[:,Xindxs] + self.vel[:,Yindxs]*self.vel[:,Yindxs])

        ## individual relative velocity
        self.relsp  = sqrt(self.u[:,Xindxs]*self.u[:,Xindxs] + self.u[:,Yindxs]*self.u[:,Yindxs])

        ## Individual radius againt to the center of mass
        self.rad   = sqrt(self.r[:,Xindxs] * self.r[:,Xindxs] + self.r[:,Yindxs] * self.r[:,Yindxs])
        self.aveR = mean(self.rad,axis=1)
        meanrad = tile(self.aveR,(self.N,1)).transpose()
        self.relrad= self.rad/meanrad

        ## Individual angle to the X axis
        self.theta = arctan2(self.r[:,Xindxs],self.r[:,Yindxs])

        ## virial of the system
        self.virial = sum( self.r[:,Xindxs]*self.vel[:,Xindxs] + self.r[:,Yindxs]*self.vel[:,Yindxs], axis=1 )

        ## normalized 
        self.Lnormed = self.angmom/(self.speed*sqrt(self.r[:,Xindxs]**2+self.r[:,Yindxs]**2))
        
    def Average_Coeff_Relative(self):
        ps = [ self.Coeff_RelativeDiffRot(i,j) for i in range(self.N) for j in range(self.N) if i != j ]
        
        return mean(ps,axis=0)

    def Coeff_RelativeDiffRot(self,i,j):
        if i == j:
            print("Error : i and j must be different")
            exit(0)

        init_dist = sqrt( (self.r[0,2*i]-self.r[0,2*j])**2 + (self.r[0,2*i+1]-self.r[0,2*j+1])**2 )
        ref_dist  = add.accumulate( tile(init_dist,self.timelen) )
        dist = sqrt( (self.r[:,2*i]-self.r[:,2*j])**2 + (self.r[:,2*i+1]-self.r[:,2*j+1])**2 )
        rot  = arctan2(self.r[:,2*i]-self.r[:,2*j], self.r[:,2*i+1]-self.r[:,2*j+1])
        diff_ij       = add.accumulate(dist) - ref_dist
        diff_ij_abs   = abs( add.accumulate(dist) - ref_dist )
        diff_ij_ratio = add.accumulate(dist)/ref_dist
        rot_ij        = add.accumulate(rot)
        #################################################################
        # self.dist_rel4init = dist/init_dist
        
        # ftdata=lf.LeastSquareOptimize(xdata, self.dist_rel4init, 'linear')
        # ftdata.optimize(1.0,0.0)
        # p_d = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        #################################################################

        xdata= arange(0,len(dist))

        ftdata=lf.LeastSquareOptimize(xdata, diff_ij, 'linear')
        ftdata.optimize(1.0,0.0)
        p = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        ftdata=lf.LeastSquareOptimize(xdata, diff_ij_abs, 'linear')
        ftdata.optimize(1.0,0.0)
        p_abs = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        ftdata=lf.LeastSquareOptimize(xdata, diff_ij_ratio, 'linear')
        ftdata.optimize(1.0,0.0)
        p_ratio = ftdata.GetOptParams()
        # ftdata.ShowFittingFunc('linear')

        ftdata=lf.LeastSquareOptimize(xdata, rot_ij, 'linear')
        ftdata.optimize(1.0,0.0)
        p_rot = ftdata.GetOptParams()

        # self.rel_diff_coeff = [p[0], p_abs[0], p_ratio[0], p_rot[0]]
        return p[0], p_abs[0], p_ratio[0], p_rot[0]

        # ts = fac.TimeSeries(dist,1)
        # ts.ShowPowerSpectorum()

        # plt.close()
        # [hist,bins] = histogram(dist,bins=100)
        # plt.bar(bins[0:len(hist)],hist,width=bins[1]-bins[0])
        # plt.show()
