##load in pretrack file, track with trackpy and calculate displacement matrix
##inputs must be of form 'filename', need single quotes

pretrackbig=input('pretrack file for big particles:') #load in pretrack file, gdf
pretracksmall=input('pretrack file for small particles:')
trackbig=input('filename to save big particle track file(.h5):') #name of file for track file to be saved as .csv
tracksmall=input('filename to save small particle track file(.h5):')
prefbig=input('prefix for position table and displacement matrix for big particles:') #prefix for position table and displacement matrix files
prefsmall=input('prefix for position table and displacement matrix for small particles:')
lifetime=input('lifetime=')

import softpy
import numpy as np
import pandas as pd
import trackpy as tp

ptbig=softpy.sftio.read_gdf(pretrackbig)
ptbig=ptbig.T
ptsmall=softpy.sftio.read_gdf(pretracksmall)
ptsmall=ptsmall.T

fbig=pd.DataFrame(ptbig,columns=['x','y','mass','size','ecc','frame']) #convert to data frame
fsmall=pd.DataFrame(ptsmall,columns=['x','y','mass','size','ecc','frame']) 

(frame for fnum, frame in fbig.groupby('frame'))
trabig=pd.concat(tp.link_df_iter((frame for fnum, frame in fbig.groupby('frame')),5,memory=10))
(frame for fnum, frame in fsmall.groupby('frame'))
trasmall=pd.concat(tp.link_df_iter((frame for fnum, frame in fsmall.groupby('frame')),5,memory=10))

trabig.to_hdf(trackbig,'tablebig') #save track file
trasmall.to_hdf(tracksmall,'tablesmall')

disbig,tabbig=softpy.sftio.tra_df_2_disp(trabig,lifetime,prefix=prefbig)
dissmall,tabsmall=softpy.sftio.tra_df_2_disp(trasmall,lifetime,prefix=prefsmall)

dismatbig=np.concatenate((disbig[:,0,:],disbig[:,1,:]),axis=1)
dismatsmall=np.concatenate((dissmall[:,0,:],dissmall[:,1,:]),axis=1)

np.savetxt(prefbig+'_displacements.txt',dismatbig) #save displacement matrix
np.savetxt(prefbig+'_positiontable.txt',tabbig) #save average position table 
np.savetxt(prefsmall+'_displacements.txt',dismatsmall) #save displacement matrix
np.savetxt(prefsmall+'_positiontable.txt',tabsmall) #save average position table 







