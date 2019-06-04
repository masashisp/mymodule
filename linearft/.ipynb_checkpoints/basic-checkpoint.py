# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from numpy import array, zeros, sqrt, pi, exp, arange, mean, double, log, cumsum, sin, divide
import numpy
from scipy.optimize import leastsq
from scipy.stats import distributions, poisson, lognorm
import scipy.special as spsp
import scipy

def histinfo(xbins,nbin):
    dx = xbins[1] - xbins[0]
    offset = 0.5 * dx
    xdatafit = xbins[0:nbin] + offset - xbins[0]
    xdataplot = xbins[0:nbin] + offset
    return dx, offset, xdatafit, xdataplot

##### Fitting distribution with the parameters without dividing
def gauss(x, p):
    """
    p[0] : mu, location
    p[1] : inverse of sigma, shape
    """
    return ( p[1]/sqrt(2.0*pi) ) * exp( -( (x-p[0])**2 ) * ( 0.5 * (p[1]**2) ) )

def dualgauss(x, p):
    """
    p[0] : mu
    p[1] : inverse of sigma
    """
    return ( ( p[1] * sqrt(2.0)*exp(-0.5*(p[1]**2)*(x[0]-p[0])**2) ) / (4.0*sqrt(pi)) ) * np.erf(x[1]-x[0] + p[0])

def twosidegauss(x, p):
    """
    p[0] : mu1
    p[1] : inverse of sigma1
    p[2] : mu2
    p[3] : inverse of sigma2
    """
    return 0.5*( p[1] * sqrt(2.0)*exp(-0.5*(p[1]**2)*(x-p[0])**2) ) + 0.5*( p[3] * sqrt(2.0)*exp(-0.5*(p[3]**2)*(x-p[2])**2) )

def gauss_cum(x, p):
    """
    p[0] : mu, location
    p[1] : inverse of sigma, shape
    """
    # xs=arange(-
    return ( p[1]/sqrt(2.0*pi) ) * exp( -( (x-p[0])**2 ) * ( 0.5 * (p[1]**2) ) )

def gammadist(x, p):
    """
    p[0] : k, shape
    p[1] : theta, scale
    # p[2] : m, location
    """
    # return ( (p[1]**p[0]) * ( x**(p[0]-1) ) * exp( -p[1]*x ) ) / spsp.gamma(p[0])
    return ( ( x**(p[0]-1) ) * exp( -x/p[1] ) ) / ( spsp.gamma(p[0]) * (p[1]**p[0]) )
    # return (((x-p[2])**(p[0]-1))*exp(-(x-p[2])/p[1]))/(spsp.gamma(p[0])*(p[1]**p[0]))

def gammafixk(x, p, k=2.0):
    """
    k : k, shape is fixed
    p[0] : theta, scale
    """
    # return ( (p[1]**p[0]) * ( x**(p[0]-1) ) * exp( -p[1]*x ) ) / spsp.gamma(p[0])
    return ( ( x**(k-1) ) * exp( -x/p[0] ) ) / ( spsp.gamma(k) * (p[0]**k) )
    # return (((x-p[2])**(p[0]-1))*exp(-(x-p[2])/p[1]))/(spsp.gamma(p[0])*(p[1]**p[0]))

def gengamma(x, p):
    """
    p[0] : a, scale
    p[1] : d, power
    p[2] : p, weibull scale

    p[1] = p[2] + 1 : weibull distribution
    p[2] = 1 : gamma distribution
    """
    # return  ((p[2]/p[0])*((x/p[0])**(p[1]-1))*exp(-(x/p[0])**p[2]))/(spsp.gamma(p[1]/p[2]))
    return  (p[2]*(x**(p[1]-1))*exp(-(x/p[0])**p[2]))/(spsp.gamma(p[1]/p[2])*p[0]**p[1])

def weibull(x, p):
    """
    p[0] : k, shape
    p[1] : inverse of nu, scale
    """
    if isinstance(x, (int, float, complex)) and x<0.0:
        ans = 0.0
    else:
        x = array(x)
        x[x<0.0] = 0.0
    if p[0] < 0.0:
        p[0] = 0.0
    ans =  (p[0]*p[1])*( (x*p[1])**(p[0]-1) )*exp( -( (x*p[1])**p[0] ) )
    return ans
    # return ( (p[0]/p[1])*(x/p[1])**(p[0]-1) )*exp(-(x/p[1])**p[0])

def maxwell(x,p):
    """
    p[0] : alpha, shape corresponding to m/2kbT
    p[1] : m, location
    """
    # if isinstance(x, (int, long, float, complex)) and x<p[2]:
    #     ans = 0.0
    # else:
    #     x = array(x)
    #     x[x<p[1]] = p[1]
    # return 4.0*pi*((p[0]/pi)**1.5)*((x-p[1])**2)*exp(-p[0]*(x-p[1])**2)  ## 3 dimension !!
    # return 2.0*p[0]*(x-p[1])*exp(-p[0]*(x-p[1])**2) ## 2 dimension !!
    if p[0] < 0.0:
        p[0] = 0.0
    return 2.0*p[0] * x * exp(-p[0]*(x**2)) ## 2 dimension !!

def MittagLeffler(data,p,itr=171):
    """
    p[0] : alpha, shape 0 < alpha < 1
    p[1] : inverse of c, scale
    itr  : iteration to calculate : spsp.gamma has 'inf' larger than 171 
    """
    data = double(data)
    if p[0] < 0.0:
        p[0] = 0.0
    elif p[0] > 1.0:
        p[0] = 1.0
    p = double(p)
    ks = double(list(range(1,itr)))
    s = zeros(len(data),dtype=double)
    # s = double(s)
    # coeff = [ ( spsp.gamma(p[0]*k) * sin(pi*p[0]*k) )/spsp.gamma(k) for k in ks ]
    # for i,x in enumerate(data):
    #     numpy.cumsum([ (-x*p[1])**(k-1)*coeff[j] for j,k in enumerate(ks) ],dtype=double)
    for i,x in enumerate(data):
        for k in ks:
            tmp = ( ((-x*p[1])**(k-1) ) * spsp.gamma(p[0]*k) * sin(pi*p[0]*k) )/spsp.gamma(k)
            # print tmp
            if tmp != tmp: # if the tmp has nan, do not add to s[i]
                break
            elif tmp == float("inf") or tmp == float("-inf"):
                break
            s[i] += double(tmp)
    return (p[1]/pi)*s

# def TwoSideML(data, p, itr=171):
#     '''
#     Two-side Mittag-Leffler Distribution
    
#     '''

def poissondist(x,p):
    return poisson.pdf(x,p)

def poisson_cum(x,p):
    return poisson.cdf(x,p)

def modifiedgamma(x,p):
    """
    p[0] : alpha
    p[1] : beta
    """
    return ( p[1] * ( x**p[0] ) * exp( -x**p[1] ) ) / spsp.gamma((p[0]+1)/p[1])

def modifiedgamma2(x,p):
    """
    p[0] : d
    p[1] : p
    p[2] : 1/a
    """
    return ( (p[1]*p[2]**p[0]) * ( x**(p[0]-1) ) * exp( -(p[2]*x)**p[1] ) ) / spsp.gamma((p[0])/p[1])

def logNorm(x, p):
    return lognorm.pdf(x, p[0], loc=p[1], scale=p[2])

def logNorm_cum(x, p):
    return lognorm.cdf(x, p[0], loc=p[1], scale=p[2])

########################################################################################
########################################################################################

def KLdist(p,q,delta):
    """
    Caclutate Kullback-Leibler distance
    p : objective probability
    q : measured probability
    """
    # p = divide(p,lengthOfData)
    # q = divide(q,lengthOfData)
    # KL = -1.0 * numpy.nansum(p * log(p/q)) #+ np.log(delta)
    KL = -1.0 * numpy.nansum(spsp.xlogy(p,p/q)) #+ log(delta)
    return KL

def SetEstimateDistrib(disttype):
    if disttype == 'gauss':
        paranames = [r'$\mu$',r'$\sigma$']
        dist = gauss
    elif disttype == 'dualgauss':
        paranames = [r'$\mu$',r'$\sigma$']
        dist = dualgauss
    elif disttype == 'twosidegauss':
        paranames = [r'$\mu_1$', r'$\sigma_1$', r'$\mu_2$', r'$\sigma_2$']
        dist = twosidegauss
    elif disttype == 'gauss_cum':
        paranames = [r'$\mu$',r'$\sigma$']
        dist = gauss_cum
    elif disttype == 'gamma':
        paranames = [r'$k$',r'$\theta$']
        dist = gammadist
    elif disttype == 'gammafixk':
        paranames = [r'$\theta$']
        dist = gammafixk
    elif disttype == 'gengamma':
        paranames = [r'$a$',r'$d$',r'$p$']
        dist = gengamma
    elif disttype == 'weibull':
        paranames = [r'$\alpha$',r'$\gamma$']
        dist = weibull
    elif disttype == 'maxwell':
        paranames = [r'$\alpha$',r'$m$']
        dist = maxwell
    elif disttype == 'mittagleffler' or disttype == 'MittagLeffler':
        paranames = [r'$\alpha$',r'$c$']
        dist = MittagLeffler
    elif disttype == 'poisson':
        paranames = [r'$¥lambda$']
        dist = poissondist
    elif disttype == 'poissonCDF':
        paranames = [r'$¥lambda$']
        dist = poisson_cum
    elif disttype == 'modifiedgamma':
        paranames = [r'$\alpha$',r'$\beta$']
        dist = modifiedgamma
    elif disttype == 'modifiedgamma2':
        paranames = [r'$d$',r'$p$',r'$a$']
        dist = modifiedgamma2
    elif disttype == 'lognorm':
        paranames = [r'$\mu$',r'$\sigma$']
        dist = logNorm
    elif disttype == 'lognormCDF':
        paranames = [r'$\mu$',r'$\sigma$']
        dist = logNorm_cum

    return dist, paranames

def GetInitParams(estmethod,disttype,data):
    if disttype == 'gauss':
        if estmethod =='leastsq':
            distpara = array([mean(data),1])
        elif estmethod =='MLE':
            distpara = None
    elif disttype == 'twosidegauss':
        if estmethod == 'leastsq':
            distpara = array([mean(data)+1.0, 1, mean(data)-1, 1])
        elif estmethod == 'MEL':
            distpara = None
    elif disttype == 'gamma':
        if estmethod =='leastsq':
            distpara = array([2,2.0]) #,mean(data)])
            # distpara = array([2,0.5]) #,mean(data)])
        elif estmethod =='MLE':
            distpara = array([2])
    elif disttype == 'gammafixk':
        if estmethod =='leastsq':
            distpara = array([2.0]) #,mean(data)])
            # distpara = array([2,0.5]) #,mean(data)])
        elif estmethod =='MLE':
            distpara = array([2])
    elif disttype == 'gengamma':
        if estmethod =='leastsq':
            distpara = array([2.,6.,6.])
        elif estmethod =='MLE':
            distpara = array([2.,6.,6.])
    elif disttype == 'weibull':
        if estmethod =='leastsq':
            if any(max(data) == data[0:10]):
                k = 1.0
            else:
                k=1.1
            distpara = array([k,1.0])
        elif estmethod =='MLE':
            distpara = array([2])
    elif disttype == 'maxwell':
        if estmethod =='leastsq':
            distpara = array([1.0, 0.0])
            # distpara = array([1.0, 1.0])
        elif estmethod =='MLE':
            distpara = None
    elif disttype == 'mittagleffler':
        if estmethod =='leastsq':
            distpara = array([0.4,1.0])
        elif estmethod =='MLE':
            distpara = array([2])
    elif disttype == 'modifiedgamma':
        if estmethod =='leastsq':
            distpara = array([2,2.0]) #,mean(data)])
        elif estmethod =='MLE':
            distpara = array([2])
    elif disttype == 'modifiedgamma2':
        if estmethod =='leastsq':
            distpara = array([2,2.0,1.0]) #,mean(data)])
        elif estmethod =='MLE':
            distpara = array([2])
    elif disttype == 'dualgauss':
        if estmethod =='leastsq':
            distpara = array([mean(abs(data)),1])
        elif estmethod =='MLE':
            distpara = None

    return distpara

class MLE(object):
    def __init__(self, data, optdist):
        """
        data    : data you want to fit, NOT FREQUENCY DATA
        optdist : a function you will fit to the data
        """
        self.data = data
        self.SetDistribution(optdist)

    def Estimate(self,params=None):
        if self.myrv.numargs==0:
            self.params = self.myrv.fit(self.data)
        elif self.myrv.numargs==1:
            self.params = self.myrv.fit(self.data,params[0],loc=params[2],shape=params[1])
        elif self.myrv.numargs==2:
            self.params = self.myrv.fit(self.data,params[0],params[1])

    def SetDistribution(self, optdist):
        if optdist=='gamma':
            self.myrv = scipy.stats.gamma
        elif optdist=='gauss':
            self.myrv = scipy.stats.norm
        elif optdist=='weibull':
            self.myrv = scipy.stats.weibull_min
        else:
            print("Choose distribution from gamma, gauss, weibull")

    def ShowFittingFunc(self, plotstyle):
        [hist, xbins] = scipy.histogram(self.data,bins=100,density=True)
        if len(self.params) == 2:
            y  = self.myrv.pdf(xbins[1:],loc=self.params[0],scale=self.params[1])
        elif len(self.params) == 3:
            y  = self.myrv.pdf(xbins[1:],self.params[0],loc=self.params[1],scale=self.params[2])
        elif len(self.params) == 4:
            y  = self.myrv.pdf(xbins[1:],self.params[0],self.params[1],loc=self.params[2],scale=self.params[3])

        self.xarray = xbins[0:100]
        self.ft_xarray = xbins[0:100]
        w = xbins[1]-xbins[0]

        if plotstyle=='linear':
            plt.bar(self.xarray, hist, width=w,label='original')
            plt.plot(self.ft_xarray, y)
            # plt.plot(self.ft_xarray, y,label=self.label)
        # elif plotstyle=='log':
        #     plt.bar(self.xarray, hist,width=w,label='original')
        #     plt.loglog(self.ft_xarray, y)
        #     # plt.loglog(self.ft_xarray, y,label=self.label)
        elif plotstyle=='semilogy':
            plt.bar(self.xarray, hist,width=w,log=True,label='original')
            plt.semilogy(self.ft_xarray, y)
            # plt.semilogy(self.ft_xarray, y,label=self.label)
        # elif plotstyle=='semilogx':
        #     plt.bar(self.xarray, hist,width=w,label='original')
        #     plt.semilogx(self.ft_xarray, y)
        #     # plt.semilogx(self.ft_xarray, y,label=self.label)
        else:
            print("Choose plot style from : linear, log, semilogy, semilogx")

        plt.legend(loc='best')
        plt.show()

class LeastSqOptForDist(object):
    def __init__(self, xarray, data, optdist, optpara = None):
        """
        xarray  : x-axis data
        data    : data you want to fit
        optfunc : a function you will fit to the data
        optpara : initialized parameter for the fitting function
        """
        self.mask = False
        self.data = array(data)
        self.xarray = xarray
        self.modeldist = optdist
        if type(optpara) == type(None):
            self.dist, self.modelpara = SetEstimateDistrib(self.modeldist.__name__)
        else:
            self.dist, _ = SetEstimateDistrib(self.modeldist.__name__)
            self.modelpara   = optpara
        self.paranum = len(self.modelpara)
        self.param_output = zeros([self.paranum, self.paranum])

    def residueWP(self, p, x, y):
        res = (y - self.dist(x, p))
        return res

    def WeibullPlotOptimize(self, params):
        self.dx=self.xarray[1]-self.xarray[0]
        self.wpxdata = log(self.xarray)
        self.wpydata = log(-log(1.0-cumsum(self.data)*self.dx))
        self.param_outputWP = leastsq(self.residueWP,
                                    params,
                                    args=(self.wpxdata,
                                          self.wpydata),
                                    full_output=True)

        self.wpk = self.param_outputWP[0][0]
        self.wpl = exp(self.param_outputWP[0][1]/self.wpk)

    def GetWeibullPlotParameters(self):
        return self.param_outputWP[0]

    def residue(self, p, x, y):
        res = (y-self.modeldist(x,p))
        return res

    def optimize_mask_data(self, init,end):
        self.m_data = self.data[init:end]
        self.ft_xarray = self.xarray[init:end]
        self.mask = True

        self.param_output = leastsq(self.residue,
                                    self.modelpara,
                                    args=(self.ft_xarray,
                                          self.m_data),
                                    full_output=True)

    def optimize(self):
        self.ft_xarray = self.xarray

        # tmp = leastsq(self.residue,
        #               self.modelpara,
        #               args=(self.ft_xarray,
        #                     self.data),
        #               full_output=True)
        # print tmp
        self.param_output = leastsq(self.residue,
                                    self.modelpara,
                                    args=(self.ft_xarray,
                                          self.data),
                                    full_output=True)
    def GetOptParams(self,disttype):
        # if disttype == 'gauss':
        #     para = [self.param_output[0][0],1.0/self.param_output[0][1]]
        # elif disttype == 'gamma':
        #     para = self.param_output[0]
        # elif disttype == 'weibull':
        #     para = [self.param_output[0][0],1.0/self.param_output[0][1]]
        # elif disttype == 'maxwell':
        #     para = self.param_output[0]
        # elif disttype == 'mittagleffler':
        #     para = self.param_output[0]
        # return para
        return self.param_output[0]

    def SetDistribution(self,dist):
        self.modeldist = dist

    def ShowFittingFunc(self, plotstyle):
        y  = self.modeldist(self.ft_xarray, self.param_output[0])

        if plotstyle=='linear':
            plt.plot(self.xarray, self.data,'.',label='original')
            plt.plot(self.ft_xarray, y)
            # plt.plot(self.ft_xarray, y,label=self.label)
        elif plotstyle=='log':
            plt.loglog(self.xarray, self.data,'.',label='original')
            plt.loglog(self.ft_xarray, y)
            # plt.loglog(self.ft_xarray, y,label=self.label)
        elif plotstyle=='semilogy':
            plt.semilogy(self.xarray, self.data,'.',label='original')
            plt.semilogy(self.ft_xarray, y)
            # plt.semilogy(self.ft_xarray, y,label=self.label)
        elif plotstyle=='semilogx':
            plt.semilogx(self.xarray, self.data,'.',label='original')
            plt.semilogx(self.ft_xarray, y)
            # plt.semilogx(self.ft_xarray, y,label=self.label)
        else:
            print("Choose plot style from : linear, log, semilogy, semilogx")

        plt.legend(loc='best')
        plt.show()

class LeastSquareOptimize(object):
    def __init__(self, xarray, data, optfunc):
        """
        xarray  : x-axis data
        data    : data you want to fit
        optfunc : Choose 'linear', 'power', 'exp', or 'gauss' to set function you will fit
        """
        self.mask = False
        self.data = array(data)
        self.xarray = xarray
        self.param_output = [[0,0],[0,0]]

        if optfunc == "linear":
            self.modelfunc = self.linearfunc
        elif optfunc == "power":
            self.modelfunc = self.powerfunc
        elif optfunc == "exp":
            self.modelfunc = self.expfunc
        elif optfunc == "gauss":
            self.modelfunc = self.gaussfunc
        elif optfunc == None:
            print("Choose function to fit. Input 'linear', 'power', 'exp', or 'gauss'")
        else:
            self.modelfunc = optfunc

    def linearfunc(self, x, p):
        self.label = f'y={self.param_output[0][0]:.3f}x+{self.param_output[0][1]:.3f}'
        func = p[0] * x + p[1]
        return func

    def powerfunc(self, x, p):
        self.label = f'y={self.param_output[0][0]:.3f}x^{self.param_output[0][1]:.3f}'
        func = p[0] * x ** p[1]
        return func

    def expfunc(self, x, p):
        self.label = f'y={self.param_output[0][0]:.3f}exp({self.param_output[0][1]:.3f}x)'
        func = p[0]*exp(p[1]*x)
        return func

    def gaussfunc(self, x, p):
        self.label = f'y={self.param_output[0][0]:.3f}exp(-x^{self.param_output[0][1]:.3f})'
        func = p[0]*exp(-x**p[1])
        return func

    def residue(self, p, x, y):
        res = (y-self.modelfunc(x,p))
        return res

    def optimize_mask_data(self,a,b,init,end):
        self.m_data = self.data[init:end]
        self.ft_xarray = self.xarray[init:end]
        self.mask = True
        p0 = [a, b]
        self.param_output = leastsq(self.residue,
                                    p0,
                                    args=(self.ft_xarray,
                                          self.m_data),
                                    full_output=True)

    def optimize(self, p):
        self.ft_xarray = self.xarray
        self.param_output = leastsq(self.residue,
                                    p,
                                    args=(self.xarray,
                                          self.data),
                                    full_output=True)

    def GetOptParams(self):
        return self.param_output[0]

    def GetAIC(self, p, mask=False):
        if mask:
            self.m_data = self.data[mask[0]:mask[1]]
            self.ft_xarray = self.xarray[mask[0]:mask[1]]
            return -2.0 * log( sum( self.residue(p, self.ft_xarray, self.m_data)**2 ) ) + 2.0 * float(len(p))
        else:
            return -2.0 * log( sum( self.residue(p, self.xarray, self.data)**2 ) ) + 2.0 * float(len(p))

    def GetBIC(self, p, mask=False):
        if mask:
            self.m_data = self.data[mask[0]:mask[1]]
            self.ft_xarray = self.xarray[mask[0]:mask[1]]
            return -2.0 * log( sum( self.residue(p, self.ft_xarray, self.m_data)**2 ) ) + log(len(self.ft_xarray)) * float(len(p))
        else:
            return -2.0 * log( sum( self.residue(p, self.xarray, self.data)**2 ) ) + 2.0 * float(len(p))

    def ShowFittingFunc(self, plotstyle):
        y  = self.modelfunc(self.ft_xarray, self.param_output[0])

        if plotstyle=='linear':
            plotter = plt.plot
        elif plotstyle=='log':
            plotter = plt.loglog
        elif plotstyle=='semilogy':
            plotter = plt.semilogy
        elif plotstyle=='semilogx':
            plotter = plt.semilogx
        else:
            print("Choose plot style from : linear, log, semilogy, semilogx")

        plotter(self.xarray, self.data,'.',label='original')
        plotter(self.ft_xarray, y, label=self.label)
        plt.legend(loc='best', frameon=False)
        plt.show()

if __name__=='__main__':
    import numpy as np

    Px = array( [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5] )
    Py = array( Px + np.tile(2.0,len(Px)) + np.random.random(len(Px)) )

    lso = LeastSquareOptimize(Px,Py,'linear')
    lso.optimize(1,2)
    print(lso.GetOptParams())
    lso.ShowFittingFunc('linear')

    Px = array( np.arange(0,10) )
    Py = array( Px**3 + np.random.random(len(Px)) )

    lso = LeastSquareOptimize(Px,Py,'power')
    lso.optimize(1,2)
    print(lso.GetOptParams())
    lso.ShowFittingFunc('log')

    # ########################################################
    # gausspara = array([10,1.5])
    # dist = gauss

    # gammapara = array([100.0,0.1])

    # weibullpara = array([2.0,0.8])

    # poissonpara = [12]
    
    # dist = gammadist
    # distpara = gammapara
    # plotstyle = 'linear'
    
    # #######################################################
    # Px = arange(0.0,20.0,0.01)
    # Py = array( dist(Px-min(Px),distpara) + np.random.random(len(Px))*0.001 )
    # [ plt.plot(Px,gammadist(Px,[15+10*i,0.1])) for i in range(10) ]
    # plt.show()
             
    # fitpara = array([5,1])
    # lsd = LeastSqOptForDist(Px,Py,dist,distpara)
    # lsd.optimize()
    # print lsd.GetOptParams()
    # lsd.ShowFittingFunc(plotstyle)

    # data = scipy.randn(10000)-3.0
    # data = scipy.stats.weibull_min.rvs(2,size=10000,loc=4,scale=0.5)
    # histdata,bins = scipy.histogram(data,bins=100,density=True)
    # weibullpara=array([2.0,0.5,1.0])

    # mle = MLE(data,'weibull')
    # mle.Estimate(weibullpara)
    # # print mle.params[0],mle.params[2],mle.params[1]
    # mle.ShowFittingFunc('linear')

    # print weibull([0.1,1,4],weibullpara)

    # lsd = LeastSqOptForDist(bins[0:100],histdata,weibull,weibullpara)
    # lsd.optimize()
    # # print lsd.GetOptParams()
    # lsd.ShowFittingFunc('linear')

    # dx = 0.02
    # # x = np.arange(0,5,dx)
    # # p1 = scipy.stats.maxwell.pdf(x)
    # # p2 = scipy.stats.norm.pdf(x,2.5,1)
    # x = np.arange(-10,10,dx)
    # mu1 = 0.0
    # s1  = 2.0
    # mu2 = 0.0
    # s2  = 1.0
    # p1 = scipy.stats.norm.pdf(x,mu1,s1)
    # p2 = scipy.stats.norm.pdf(x,mu2,s2)

    # ns = 10000000
    # mu = mu1
    # sigma = s1
    # y = mu + sigma*np.random.randn(ns)
    # h = np.histogram(y,bins=len(x),range=(x[0],x[-1]),density=True)
    # q = h[1][1:]-h[1][:-1]
    # p = h[0] # /(float(len(y))*q[0]) <-- equivarent to take density=True
    # lsd = LeastSqOptForDist(h[1][0:-1],p,gauss,[mu,sigma])
    # lsd.optimize()
    # params = lsd.GetOptParams()
    # print params
    # p3 = scipy.stats.norm.pdf(h[1][0:-1],loc=params[0],scale=params[1])
    # # print sum(p*q),q[0]

    # print sum(p1*dx),sum(p2*dx),sum(p3*dx),sum(p)
    # KLmy = KLdist(p1*dx,p2*dx,dx)
    # KLmy2 = KLdist(p*q[0],p2*dx,q[0])
    # KLmy3 = KLdist(p3*dx,p2*dx,q[0])
    # KLscipy = scipy.stats.entropy(p1,p2)
    # # print scipy.stats.entropy(p1,p2), scipy.stats.entropy(p1*dx,p2*dx)
    # KLscipy2 = scipy.stats.entropy(p,p2)
    # KLscipy3 = scipy.stats.entropy(p3,p2)

    # plt.plot(x,p1,label=r'scipy $\mu=$'+str(mu1)+r', $\sigma=$'+str(s1))
    # plt.plot(x,p2,label=r'scipy $\mu=$'+str(mu1)+r', $\sigma=$'+str(s2))
    # plt.plot(h[1][1:],p,label='emprical')
    # plt.plot(h[1][1:],p3,label='fitting')
    # plt.legend(loc='best')
    # plt.show()

    # print KLmy, KLmy2, KLmy3, KLscipy, KLscipy2, KLscipy3
