# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 16:30:53 2015

Inputs: eigenvectors eigenvec
@author: mgratale
"""

import numpy as np

def participation_ratio(eigenvec):
    """
    Calculate participation ratio from eigenvectors of DOS calculation
    
    inputs:
    eigenvec = eigenvectors of system

    """
    N=len(eigenvec[:,0])/2.0
    ev2=eigenvec*eigenvec  #numerator
    ev4=ev2*ev2  #denominator
    
    numer=np.sum(ev2,axis=0)**2
    denom=N*np.sum(ev4,axis=0)
    
    partrat=numer/denom
    
    return partrat
    
    
    