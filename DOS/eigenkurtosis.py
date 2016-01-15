# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:04:08 2015

@author: mgratale
"""

import numpy as np
import scipy.stats as sps

def eigenkurtosis(disp,eigvec,eigenval):
    """
    calculate kurtosis of eigenmodes
    inputs: displacement matrix/array disp (T x dim x N)
            eigvec=eigenvectors
            eigenval=eigenvalues
    """        
    [numt,dim,numpart]=disp.shape
    nfreq=len(eigenval)
    
    dxdy=np.concatenate((disp[:,0,:],disp[:,1,:]),axis=1)        
    
    P=dxdy.dot(eigvec)
    kurt=sps.kurtosis(P,axis=0)
        
    return kurt