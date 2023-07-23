# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

import numpy as np
import scipy.linalg as la
from SynchPhaseTimes2 import SynchPhaseTimes2

def PiCA(x,peaks,varargin = []):

    if varargin:
        peaks2 = varargin[1]
        flag = 1
    else:
        flag = 0

    [T0,T1] = SynchPhaseTimes2(peaks)

    A = np.dot(x[:,T0],np.transpose(x[:,T1]))
    B = np.dot(x[:,T0],np.transpose(x[:,T0]))
    A = (A + np.transpose(A)) / 2
    B = (B + np.transpose(B)) / 2

    if flag == 0:
        [V,D] = la.eigh(A,B) # this won't work
    elif flag == 1:
        [T0,T1] = SynchPhaseTimes2(peaks2)
        AA = np.dot(x[:,T0],np.transpose(x[:,T1]))
        AA = (AA + np.transpose(AA)) / 2
        [V,D] = la.eigh(A - AA,B) # eig doesn't give same result

    print(A)
    print(B)
    
    d = np.diag(D)
    I = np.argsort(d)
    I = I[::-1]
    
    print(V)
    print(I)
    print(V[:,I])

    W = np.transpose(V[:,I])
    A = np.linalg.pinv(W)

    y = np.dot(W,x)

    return(y,W,A)