# The .py implementation of TikhonovRegularization.mat (from OSET/Tools/Generic) that matches numerically the reconstruction.
# For very long signals, the processing speed can be increased by decression the lsqr tolerance. In the current setting it is chosen to match Matlab.

import numpy as np
from scipy import sparse
from scipy import linalg
from scipy.sparse.linalg import lsqr
def TikhonovRegularizationMD(x, DiffOrderOrFilterCoefs, Lambda):
    x = np.array(x)
    if(np.isscalar(DiffOrderOrFilterCoefs) == True):
        l = np.zeros(DiffOrderOrFilterCoefs, dtype=int).tolist()
        l.append(1)
        l.extend(np.zeros(DiffOrderOrFilterCoefs, dtype=int).tolist())
        h = np.diff(l,DiffOrderOrFilterCoefs).tolist()
    else:
        h = DiffOrderOrFilterCoefs
    L = len(h)
    N = x.shape[len(x.shape)-1]

    a = np.zeros(N-L, dtype=int).tolist()
    a.insert(0, h[0])
    b = h.copy()
    b.extend(np.zeros(N-L, dtype=int).tolist())

    #from scipy import linalg
    from scipy import sparse
    from scipy.sparse.linalg import lsqr
    #from scipy.sparse.linalg import inv

    D = sparse.csc_matrix(linalg.toeplitz(a,b)) 
    F = sparse.csc_matrix(Lambda*(D.T*D) + np.eye(N))

    y = lsqr(F.T, x.T, btol = 1e-14,atol = 1e-14)[0].T
    return y