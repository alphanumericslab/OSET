# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import numpy as np


def SynchPhaseTimes2(peaks):
    I = np.where(peaks)
    I = I[0]
    D = I[1:] - I[0:-1]

    if len(I) < 3:
        T0 = []
        T1 = []
    else:
        start = I[0]
        stop = I[-2]

        T1 = np.zeros(stop - start + 1)
        k = 0

        for t in range(start, stop + 1):
            T1[t - start] = I[k + 1] + round((t - I[k]) * D[k + 1] / D[k])

            if t >= I[k + 1]:
                k = k + 1

    T1 = T1.astype(int)
    T0 = np.arange(start, stop + 1)

    return (T0, T1)


peaks = np.array([0, 0, 1, 0, 0, 0, 1, 0, 1, 0])
[T0, T1] = SynchPhaseTimes2(peaks)
