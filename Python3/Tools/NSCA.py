# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

import numpy as np
import scipy.linalg as la

def NSCA(x,I,J):

    B = np.cov(x[:,I])
    C = np.cov(x[:,J])

    C = (C + np.transpose(C)) / 2
    B = (B + np.transpose(B)) / 2

    D,V = la.eig(B,C) #eig doesn't work for getting same eigenvalues as matlab eig
    print(D)
    print(V)

    d = np.diag(D)
    II = np.argsort(d)
    II = II[::-1]

    W = np.transpose(V[:,II])
    A = np.linalg.pinv(W)

    y = np.real(np.dot(W,x))

    return(y,W,A,B,C)