# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import numpy as np


def PolyFit(x, fs, N):
    x = np.transpose(x)

    L = len(x)
    n = np.arange(0, L)
    t = n / fs

    T = np.zeros((N, L))
    TX = np.zeros((N, 1))
    for i in range(0, N):
        T[i, :] = t ** (i)
        TX[i] = np.dot(x, np.transpose(T[i, :]))

    TT = np.dot(T, np.transpose(T))

    p = np.dot(np.linalg.pinv(TT), TX)

    y = np.dot(np.transpose(p), T)

    print(y)
    print(p)

    return (y, p)
