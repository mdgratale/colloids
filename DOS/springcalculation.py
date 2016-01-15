# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 21:07:21 2015

@author: mgratale
"""

import numpy as np

def springcalculation(kmat,x,y,mass):
    """
    kmat = spring constant matrix
    x,y,mass = x positions, y positions, and masses from average position table
    """
    
    numpart=int(len(x))
    xmat=np.ones((numpart,numpart))*x
    ymat=np.ones((numpart,numpart))*y
    dx=xmat-xmat.T
    dy=ymat-ymat.T
    dr=np.sqrt((dx*dx)+(dy*dy))

    m1=np.ones((numpart,numpart))*mass
    mm2=m1*m1.T
        
    kxx=kmat[0:numpart,0:numpart]
    kyy=kmat[numpart:2*numpart,numpart:2*numpart]
    kmag=-1*(kxx+kyy)
    
    springs=np.zeros((numpart,numpart,3))
    springs[:]=np.NaN
    springs[:,:,0]=dr
    springs[:,:,1]=mm2
    springs[:,:,2]=kmag

    return springs    
    