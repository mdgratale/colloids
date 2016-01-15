"""
Function to correct for drift in displacement array
"""

import numpy as np

def mgdisp_correction(disp):

	newdx=disp.copy()

	[nt,dim,numpart]=disp.shape

	for i in range(0,nt):
         idxi=disp[i,0,:]!=0
         idyi=disp[i,1,:]!=0
         newdx[i,0,idxi]=disp[i,0,idxi]-np.mean(disp[i,0,idxi])
         newdx[i,1,idyi]=disp[i,1,idyi]-np.mean(disp[i,1,idyi])

	return newdx
