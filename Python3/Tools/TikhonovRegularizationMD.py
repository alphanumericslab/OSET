import numpy as np
#from scipy.linalg import toeplitz, 
from scipy import linalg
from scipy import sparse
def TikhonovRegularization(x, DiffOrderOrFilterCoefs, Lambda):
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
    
    D = sparse.csr_matrix(linalg.toeplitz(a,b)) 
    F = (Lambda*(D.T*D) + np.eye(N))
    
    # Both ways of performing the matrix right division are equally efficient to my knowledge.
    # If there is any reason to chose one, go for it.
    y = np.linalg.lstsq(F.T,x.T,rcond=None)[0].T
    #y = np.dot(x, linalg.pinv(F)) 
    return y