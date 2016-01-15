# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 10:08:08 2015

@author: mgratale
"""

def eigenvalue_extrapolation(eg1,eg2,eg3,eg4,t1,t2,t3,t4):
    """
    Extrapolate eigenvalues to infinite time, currently only written for 4 time frame inputs

    Inputs:  
    eg1 = first set of eigenvalues (longest time)
    eg2 = second set of eigenvalues (second longest time)
    eg3 = third set of eigenvalues (next longest time)
    eg4 = fourth set of eigenvalues (shortest time)
    t1,t2,t3,t4 = number of frames for each set of eigenvalues           

    """    
    nfreq=len(eg1) #number of frequencies/degrees of freedom
    times=np.array([t1,t2,t3,t4])
    R=float(nfreq)/times
    
    eg1=np.sort(eg1)
    eg2=np.sort(eg2)
    eg3=np.sort(eg3)
    eg4=np.sort(eg4)
    
    freq=np.sqrt(np.array([eg1,eg2,eg3,eg4]))
    newfreq=eg1.copy()
    newfreq[:]=np.NAN
    slopes=newfreq.copy()        
    
    for i in range(0,nfreq):
        m,b=np.polyfit(R,freq[:,i],1)
        newfreq[i]=b
        slopes[i]=m
        
    return newfreq,slopes
        