# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:41:41 2015

@author: mgratale
"""

import numpy as np

def bigsmall_partfrac(eigenvec,masses,thresh=1):
    """
    function to calculate participation fraction of big and small particles in colloidal glass
    eigenvec=eigenvectors
    masses=1D array of particles masses from avg. position table
    thresh=cutoff for big and small particles in masses
    """
    
    numpart=len(eigenvec[:,0])/2
    amp=eigenvec[0:numpart,:]**2 + eigenvec[numpart:2*numpart,:]**2    
    
    big=masses>thresh
    small=masses<=thresh
    
    bpf=np.sum(amp[big,:],axis=0)
    spf=np.sum(amp[small,:],axis=0)
    
    return bpf,spf
    
    