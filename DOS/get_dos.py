""""
Calculate density of states from displacement matrix
"""
import numpy as np

def get_dos(disp,masses=None,thresh=0,mass1=1.0,filename=None):
    """
     disp=(time,dim,#particles) displacement matrix
     masses=array of masses for particles, left as None for monodispersed system
     thresh=threshhold for mass ratio, mass1=mass ratio
    """

    numpart=len(disp[0,0,:])

    dmx=np.concatenate((disp[:,0,:],disp[:,1,:]),axis=1) #turn (T x dim x N) array into 2D (T x 2N) array

    if masses is None:
        masses=np.zeros(shape=(1,len(disp[0,0,:])))

    masses[masses>=thresh]=mass1
    masses[masses<thresh]=1.0
    m1=np.zeros((1,2*numpart))
    m1[0,0:numpart]=masses
    m1[0,numpart:2*numpart]=masses
    m=np.ones((2*numpart,2*numpart))*m1
    mmask=m*(m.T)
	
    print('Calculating covariance matrix')
    covmat=np.cov(dmx,rowvar=0) #covariance matrix
	
    print('Inverting to spring matrix')
    kmat=np.linalg.inv(covmat) #invert covariance matrix to find spring matrix
	
    dynmat=kmat/np.sqrt(mmask)
	
    print('Calculating eigenvalues and eigenvectors')
    w,v=np.linalg.eig(dynmat)
	
    if filename is not None:
        print('Saving data')
        np.save(filename+'_covariancematrix.npy',covmat)
        np.save(filename+'_springmatrix.npy',kmat)
        np.save(filename+'_dynamicalmatrix.npy',dynmat)
        np.save(filename+'_eigenvalues.npy',w)
        np.save(filename+'_eigenvectors.npy',v)
 
    return w,v	
		
