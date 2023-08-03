# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import numpy as np


def AdaptiveFilter(x, delay, taps, mu):
    xref = np.concatenate((np.zeros(delay + taps - 1), x[0:-delay]))
    w = np.ones(taps)

    n = np.zeros(len(x))
    e = np.zeros(len(x))

    for i in range(len(x)):
        if i > taps - 1:
            xr = xref[i : i - taps : -1]
        else:
            xr = np.concatenate((xref[i::-1], np.zeros(taps - i - 1)))

        n[i] = np.dot(w, xr)
        e[i] = x[i] - n[i]
        w += 2 * mu * e[i] * xr

    return (e, n)
