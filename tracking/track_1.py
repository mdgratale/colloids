##load in pretrack file, track with trackpy and calculate displacement matrix
##inputs must be of form 'filename', need single quotes

pretrack=input('pretrack file:') #load in pretrack file, gdf
trackfile=input('filename to save track file(.h5):') #name of file for track file to be saved as .h5
pref=input('prefix for position table and displacement matrix:') #prefix for position table and displacement matrix files
lifetime=input('lifetime=')

import softpy
import numpy as np
import pandas as pd
import trackpy as tp

pt=softpy.sftio.read_gdf(pretrack)
pt=pt.T

f=pd.DataFrame(pt,columns=['x','y','mass','size','ecc','frame']) #convert to data frame

(frame for fnum, frame in f.groupby('frame'))
tra=pd.concat(tp.link_df_iter((frame for fnum, frame in f.groupby('frame')),5,memory=10))

tra.to_hdf(trackfile,'table') #save track file

dis,tab=softpy.sftio.tra_df_2_disp(tra,lifetime,prefix=pref)

dismat=np.concatenate((dis[:,0,:],dis[:,1,:]),axis=1)

np.savetxt(pref+'displacements.txt',dismat) #save displacement matrix
np.savetxt(pref+'positiontable.txt',tab) #save average position table 


