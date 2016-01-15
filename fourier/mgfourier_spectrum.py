# -*- coding: utf-8 -*-
"""
find fourier spectrum of displacement matrix/array

@author: mgratale
"""
import numpy as np

def mgfourier_spectrum(dis):
    
    [nt,dim,partnum]=dis.shape
    dx=dis[:,0,:].copy()
    dy=dis[:,1,:].copy()
    
    fx=np.fft.fft(dx,axis=0)
    fy=np.fft.fft(dy,axis=0)
    
    f=np.concatenate((fx,fy),axis=1)
    I=np.mean(np.absolute(f),axis=1)
    
    q=np.fft.fftfreq(len(dx[:,0]))
    #q=np.linspace(1,nt,nt)/float(nt)

    return I,q    
    
    