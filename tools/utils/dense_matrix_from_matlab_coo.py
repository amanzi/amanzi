import numpy as np
import scipy.sparse as sp

def denseMatrixFromCOO(coo, shift=1):
    sparse = sp.coo_matrix( ( coo[:,2], (coo[:,0]-shift, coo[:,1]-shift)) )
    return sparse.todense()

def denseMatrixFromMatlabFile(fname, shift=1):
    coo = np.loadtxt(fname)
    return denseMatrixFromCOO(coo, shift)

