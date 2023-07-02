# Copyright (c) 2021, David Williams, Reza Sameni
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

import numpy as np
import scipy.io as spio

def KFNotch2(x,f0,fs,Qbar = [],Wlen = []):

    # input arguments
    if not Qbar:
        Qbar = .0001*max([abs(i) for i in x])

    if not Wlen:
        Wlen = int(np.ceil(.1*fs))

    # Kalman filter parameters
    w = 2 * np.pi * f0 / fs

    Wmean = 0
    Vmean = np.mean(x)
    X0 = x[1::-1]
    X0 = np.transpose([X0])
    P0 = 2*Qbar

    R = 1
    Q = Qbar
    L = 2
    A = np.array([[2*np.cos(w),-1],[1,0]])
    H = np.array([1,0])
    Xminus = X0
    Pminus = P0
    Samples = len(x)

    mem1 = np.ones(Wlen)
    mem2 = np.ones(Wlen)
    Xhat = np.zeros((2,Samples))
    innovations = np.zeros(Samples)
    Phat = np.zeros((2,2,Samples))
    Xbar = np.zeros((2,Samples))
    Pbar = np.zeros((2,2,Samples))
    Kgain = np.zeros((2,Samples))
    B = np.array([[1],[0]])

    # Forward filtering stage
    for k in range(Samples):
        Xbar[:,k] = np.transpose(Xminus)
        Pbar[:,:,k] = np.transpose(Pminus)

        Yminus = np.dot(H,Xminus) + Vmean

        K = np.divide(np.dot(Pminus,np.transpose([H])),(np.dot(np.dot(H,Pminus),np.transpose([H])) + R))
        Pplus = np.dot(np.dot((np.eye(L)-np.dot(K,[H])),Pminus),np.transpose(np.eye(L)-K*H)) + np.dot(np.dot(K,R),np.transpose(K))
        innovations[k] = x[k] - Yminus
        Xplus = Xminus + K*innovations[k]

        mem1 = np.append(np.dot(np.dot(H,Pminus),np.transpose(H)) + R,mem1[0:-1])
        mem2 = np.append(innovations[k]**2,mem2[0:-1])

        Pplus = (Pplus + np.transpose(Pplus))/2

        Xminus = np.dot(A,Xplus) + Wmean
        Pminus = np.dot(np.dot(A,Pplus),np.transpose(A)) + np.dot(B*Q,np.transpose(B))

        mu = np.mean(mem2/mem1)
        Q = Qbar * mu

        Xhat[:,k] = np.transpose(Xplus)
        Phat[:,:,k] = np.transpose(Pplus)
        Kgain[:,k] = np.reshape(K,2)

    # Backward smoothing stage
    PSmoothed = np.zeros(np.shape(Phat))
    X = np.zeros(np.shape(Xhat))
    PSmoothed[:,:,Samples-1] = Phat[:,:,Samples-1]
    X[:,Samples-1] = Xhat[:,Samples-1]

    for k in reversed(range(Samples-1)):
        S = np.dot(np.dot(Phat[:,:,k],np.transpose(A)),np.linalg.inv(Pbar[:,:,k+1]))
        X[:,k] = Xhat[:,k] + np.dot(S,(X[:,k+1] - Xbar[:,k+1]))
        PSmoothed[:,:,k] = Phat[:,:,k] - np.dot(np.dot(S,(Pbar[:,:,k+1] - PSmoothed[:,:,k+1])),np.transpose(S))

        PSmoothed[:,:,k] = (PSmoothed[:,:,k] + np.transpose(PSmoothed[:,:,k]))/2

    y1 = x - np.dot(H,Xhat)
    y2 = x - np.dot(H,X)

    return(y1,y2,Pbar,Phat,PSmoothed,Kgain)