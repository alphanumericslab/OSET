# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

import numpy as np


def RLS(x, d, p, lamb, delta):
    x = np.transpose(x)
    d = np.transpose(d)
    T = len(x)
    P = np.zeros((p, p, T))
    K = np.zeros((p, T))
    w = np.zeros((p, T))
    y = np.zeros(T)

    P[:, :, 0] = np.eye(p) / delta
    gamma = 1 / lamb

    for n in range(p - 1, T):
        if n < p:
            u = x[n::-1]
        else:
            u = x[n : n - p : -1]

        u = np.transpose(u)
        K[:, n] = (
            gamma
            * np.dot(P[:, :, n - 1], u)
            / (1 + gamma * np.dot(np.dot(np.transpose(u), P[:, :, n - 1]), u))
        )
        y[n] = d[n] - np.dot(w[:, n - 1], u)
        w[:, n] = w[:, n - 1] + K[:, n] * y[n]
        P[:, :, n] = (
            gamma * P[:, :, n - 1]
            - gamma * (K[:, n] * np.transpose(u)) * P[:, :, n - 1]
        )

    e = x - y

    return (y, e, P)
