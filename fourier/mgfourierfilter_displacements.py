# -*- coding: utf-8 -*-
"""

@author: mgratale
"""
import numpy as np

def mgfourierfilter_displacements(dis,freq,nfreq=2,delta=0.0005,filename=None):
    """
    Filter frequencies from fourier transform of displacement matrix
    Inputs:
    dis=displacement matrix/array (T x dim x N)
    freq=frequency to remove
    nfreq=number of peaks/frequencies to remove
    delta=range around peak to remove
    """
    
    [nt,dim,numpart]=dis.shape
    dx=dis[:,0,:].copy() #x displacements
    dy=dis[:,1,:].copy() #y displacements
    ref=(abs(dx)+abs(dy))>0 #reference array to find where particles don't exist       

    fx=np.fft.fft(dx,axis=0) #fft of x displacements
    fy=np.fft.fft(dy,axis=0) #fft of y displacements
    ffd=np.zeros((nt,dim,numpart),dtype=complex)
    ffd[:,0,:]=fx
    ffd[:,1,:]=fy    
    
    q=np.fft.fftfreq(len(dx[:,0])) #frequencies   
        
    k1=np.linspace(1,nfreq,nfreq)*freq
    k2=-1.0*k1
    k=np.concatenate((k1,k2),axis=1)
    
    p=np.zeros((2*nfreq,nt),dtype=bool)    
    
    for i in range(2*nfreq):
        p[i,:]=( q > k[i]-delta) & (q < delta+k[i])

    psum=np.sum(p,axis=0)
    
    good=psum==0
    good2=np.array([good])

    ffdfilter=ffd.copy()
    ffdfilter[:,0,:]=fx*good2.T
    ffdfilter[:,1,:]=fy*good2.T

    newdx=np.fft.ifft(ffdfilter[:,0,:],axis=0)
    newdy=np.fft.ifft(ffdfilter[:,1,:],axis=0)
    newdis=dis.copy()
    newdis[:,0,:]=np.real(newdx)*ref
    newdis[:,1,:]=np.real(newdy)*ref

    if filename is not None:
        np.save(filename+'displacements_fft.npy',ffd)
        np.save(filename+'displacements_fftfilter.npy',ffdfilter)
        np.save(filename+'displacements_filtered.npy',newdis)

    return newdis,ffdfilter,ffd    
        
        
    
    