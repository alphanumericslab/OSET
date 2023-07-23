# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

import numpy as np
import scipy.linalg as sp

def TikhonovRegularization(x,DiffOrderOrFilterCoefs,lamb):

    if np.isscalar(DiffOrderOrFilterCoefs):
        h = np.diff(np.concatenate((np.zeros(DiffOrderOrFilterCoefs),[1],np.zeros(DiffOrderOrFilterCoefs))),DiffOrderOrFilterCoefs)
    else:
        h = DiffOrderOrFilterCoefs

    L = len(h)
    N = np.ma.size(x,1)
    D = sp.toeplitz(np.concatenate(([h[0]],np.zeros(N-L))),np.concatenate((h,np.zeros(N-L))))
    F = np.dot(np.dot(lamb,np.transpose(D)),D) + np.eye(N)

    y = np.matmul(x, np.linalg.inv(F))

    return(y)