#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
# import numpy

from numpy import array, zeros, sqrt, pi, exp, arange, mean, double, log, cumsum, sin, divide, c_

import scipy
import scipy.special as spsp

import Pycluster as pycl
from scipy.cluster.vq import whiten

from scipy.optimize import leastsq
from scipy.stats import distributions, poisson

import sys
sys.path.append('/Users/masashi/Research/collmot/script_python')
import SwarmData as sd
from mymodule import *

def GetClustersID(knum, N, data, **kwargs):
    """
    ** argments **
    knum : number of expected clusters 
    N    : number of the particles
    data : data set, which is classified based on the clustering id of extenddata

    ** key words ** 
    extenddata : set of data that is used to calculate clustering.
                 In general, this is extended the additional features to orignal data
    weight     : define the weight of features. if weight==None, it is equals
    initit     : this is initial id to classify
    """
    extenddata = kwargs['extenddata'] if kwargs.has_key('extenddata') else data
    weight1    = kwargs['weight']     if kwargs.has_key('weight')     else None
    initid     = kwargs['initid']     if kwargs.has_key('initid')     else None

    stds = np.sqrt(np.var(extenddata,axis=0))
    whitened = whiten(extenddata)
    label,_,_=pycl.kcluster(whitened,nclusters=knum,weight=weight1,initialid=initid)

    return array(label)

def KClusters(knum, N, data, **kwargs):
# def KClusters(knum, N, extenddata, data, weight1=None, initid=None):
    """
    ** argments **
    knum : number of expected clusters 
    N    : number of the particles
    data : data set, which is classified based on the clustering id of extenddata

    ** key words ** 
    extenddata : set of data that is used to calculate clustering.
                 In general, this is extended the additional features to orignal data
    weight     : define the weight of features. if weight==None, it is equals
    initit     : this is initial id to classify
    """
    extenddata = kwargs['extenddata'] if kwargs.has_key('extenddata') else data
    weight1    = kwargs['weight']     if kwargs.has_key('weight')     else None
    initid     = kwargs['initid']     if kwargs.has_key('initid')     else None

    stds = np.sqrt(np.var(extenddata,axis=0))
    whitened = whiten(extenddata)
    label,_,_=pycl.kcluster(whitened,nclusters=knum,weight=weight1,initialid=initid)
    centroid,_=pycl.clustercentroids(whitened,clusterid=label)

    C = []
    for k in xrange(knum):
        C.append(data[label==k])

    stds = np.sqrt(np.var(extenddata,axis=0))
    centroid *= stds
    means = np.mean(data,axis=0)
    centroid[:,0] += means[0] # x
    centroid[:,1] += means[1] # y
    
    return C, centroid
    # return C1,C2,C3,C4,centroid

def sample(N,v0,Ca,Cr,la,lr,ens,init,dirname='./'):
    Xindxs = arange(0,2*N,2)
    Yindxs = arange(1,2*N,2)

    alpha = 2
    
    # fnbody = [ MakeFilenameBody(N,v0,Ca,C*Ca,la,l*la,alpha,init) for C in Cs ]

    fname = lambda Cr, e : dirname+MakeFilenameBody(N,v0,Ca,Cr,la,lr,alpha,init)+'-State-ensemble-'+str(e)+'.txt'
    sm = sd.SwarmData(fname(Cr,ens))

    knum=4
    
    for t in range(10):
        data = sm.ps_data[t,:].reshape(N,4)
        # C = sm.ps_data[t,:].reshape(N,4)
        # C = c_[sm.r[t,:].reshape(N,2),sm.angmom[t,:]]
        exdata = c_[sm.r[t,:].reshape(N,2),sm.u[t,:].reshape(N,2),sm.Lnormed[t,:]]
        # C = c_[sm.r[t,:].reshape(N,2),sm.u[t,:].reshape(N,2),sm.angmom[t,:]]
        # C = c_[sm.r[t,:].reshape(N,2),sm.u[t,:].reshape(N,2),sm.angmom[t,:]]
        # C1,C2,C3,C4 = KClusters(2, N, C)
        Cs,centroid = KClusters(knum, N, data, extenddata=exdata, weight=[1,1,1,1,2])
        # C1,C2,C3,C4,centroid = KClusters(4, N, C, sm.ps_data[t,:].reshape(N,4))

        plt.quiver(sm.pos[t,Xindxs], sm.pos[t,Yindxs],sm.vel[t,Xindxs],sm.vel[t,Yindxs])
        lt = ['ro','bo','go','yo']
        for k in xrange(knum):
            data = array(zip(*Cs[k]))
            plt.plot(data[0],data[1],lt[k])
        data = zip(*centroid)
        plt.plot(data[0],data[1],'x')

        plt.show()

def main():
    # ２種類の２次元正規分布に基づきサンプリングして観測データを生成
    mean1 = [0,0]
    cov1 = [[4,0],[0,100]]
    N1 = 1000
    X1 = np.random.multivariate_normal(mean1,cov1,N1)
    mean2 = [10,-10]
    cov2 = [[1,20],[20,50]]
    N2 = 1000
    X2 = np.random.multivariate_normal(mean2,cov2,N2)
    X = np.concatenate((X1,X2))
    # print np.shape(X)
    
    # 描画
    x,y = X.T
    plt.plot(x,y,'k.'); plt.axis('equal'); plt.show()

    print len(X)
    # kmeans2でクラスタリング
    whitened = whiten(X) # 正規化（各軸の分散を一致させる）
    centroid, label = kmeans2(whitened, k=2) # kmeans2
    C1 = []; C2 = [] # クラスタ保存用
    for i in range(len(X)):
        if label[i] == 0:
            C1 += [whitened[i]]
        elif label[i] == 1:
            C2 += [whitened[i]]
     
    # 描画
    x,y = zip(*C1)
    plt.plot(x, y, 'r.')
    x,y = zip(*C2)
    plt.plot(x, y, 'g.')
    x,y = centroid.T
    plt.plot(x, y, 'bx')
    plt.axis('equal')
    plt.show()

if __name__ == '__main__':
    # main()
    
    dirname = '../StateData/Para_C_l_0.5_N_40_t_MRK_ens_10_random_state/'
    sample(40,2,5.0,3.0,1.0,0.5,0,'random',dirname)
