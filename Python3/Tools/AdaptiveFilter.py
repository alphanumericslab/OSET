import numpy as np

def AdaptiveFilter(x,delay,taps,mu):

    xref = np.concatenate((np.zeros(delay + taps - 1),x[0:-delay]))
    w = np.ones(taps)

    n = np.zeros(len(x))
    e = np.zeros(len(x))

    for i in range(len(x)):
        if i > taps - 1:
            xr = xref[i:i-taps:-1]
        else:
            xr = np.concatenate((xref[i::-1],np.zeros(taps-i-1)))

        n[i] = np.dot(w,xr)
        e[i] = x[i] - n[i]
        w += 2 * mu * e[i] * xr
        
    return(e,n)

x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
delay = 3
taps = 5
mu = .005
[ECG_estimate, Noise_estimate] = AdaptiveFilter(x,delay,taps,mu)