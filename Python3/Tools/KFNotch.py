import numpy as np
import scipy.io as spio

def KFNotch(x,f0,fs,Q = [],R = [],gamma = 1):

    # input arguments
    if not Q:
        Q = .0001*max([abs(i) for i in x])

    if not R:
        R = np.var(x)

    # Kalman filter parameters
    w = 2 * np.pi * f0 / fs

    Wmean = 0
    Vmean = 0
    X0 = np.array(x[1::-1])
    X0 = np.transpose([X0])
    P0 = 1000*Q

    L = 2
    A = np.array([[np.cos(w),-np.sin(w)],[np.sin(w),np.cos(w)]])
    H = np.array([1,1])
    B = np.array([[1],[1]])
    BQBT = np.array([[Q,0],[0,Q]])

    Xminus = X0
    Pminus = P0
    Samples = len(x)

    VarWinlen = np.ceil(fs/10)
    mem2 = np.zeros((int(VarWinlen),1)) + R
    Xhat = np.zeros((2,Samples))
    innovations = np.zeros(Samples)
    Phat = np.zeros((2,2,Samples))
    Xbar = np.zeros((2,Samples))
    Pbar = np.zeros((2,2,Samples))
    Kgain = np.zeros((2,Samples))

    # Forward filtering stage
    for k in range(Samples):
        Xbar[:,k] = np.transpose(Xminus)
        Pbar[:,:,k] = np.transpose(Pminus)

        Yminus = np.dot(H,Xminus) + Vmean
        Yminus = Yminus[0]

        K = np.divide(np.dot(Pminus,np.transpose([H])),(np.dot(np.dot(H,Pminus),np.transpose([H])) + R))
        Pplus = np.dot(np.dot((np.eye(L)-np.dot(K,[H])),Pminus),np.transpose(np.eye(L)-K*H)) + np.dot(np.dot(K,R),np.transpose(K))
        innovations[k] = x[k] - Yminus
        Xplus = Xminus + K*innovations[k]

        mem2 = np.insert(mem2[0:-1],0,innovations[k]**2)
        mem2 = np.transpose([mem2])
        R = gamma * R + (1-gamma)*np.mean(mem2)

        Pplus = (Pplus + np.transpose(Pplus))/2

        Xminus = np.dot(A,Xplus) + Wmean
        Pminus = np.dot(np.dot(A,Pplus),np.transpose(A)) + BQBT

        Xhat[:,k] = np.transpose(Xplus)
        Phat[:,:,k] = np.transpose(Pplus)
        Kgain[:,k] = np.reshape(K,2)

    # Backward smoothing stage
    PSmoothed = np.zeros(np.shape(Phat))
    X = np.zeros(np.shape(Xhat))
    PSmoothed[:,:,Samples-1] = Phat[:,:,Samples-1]
    X[:,Samples-1] = Xhat[:,Samples-1]

    for k in reversed(range(Samples-2)):
        S = np.dot(np.dot(Phat[:,:,k],np.transpose(A)),np.linalg.inv(Pbar[:,:,k+1]))
        X[:,k] = Xhat[:,k] + np.dot(S,(X[:,k+1] - Xbar[:,k+1]))
        PSmoothed[:,:,k] = Phat[:,:,k] - np.dot(np.dot(S,(Pbar[:,:,k+1] - PSmoothed[:,:,k+1])),np.transpose(S))

        PSmoothed[:,:,k] = (PSmoothed[:,:,k] + np.transpose(PSmoothed[:,:,k]))/2

    y1 = x - np.dot(H,Xhat)
    y2 = x - np.dot(H,X)

    return(y1,y2,Pbar,Phat,PSmoothed,Kgain)