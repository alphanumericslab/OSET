import numpy as np

def filter(b,a,x):
    # creates a filter through rational transfer function with numerator b, denominator a, and input x
    
    # normalize a[0] to 1
    b = np.divide(b,a[0])
    a = np.divide(a,a[0])
    
    # create loop to run through transfer function
    y = []
    for n in range(len(x)):
        k = 0
        y = np.append(y,0)

        while k <= n:
            ind = n - k

            if ind < len(b):
                y[n] += (b[ind]*x[k])

            if ind < n and k < len(a):
                y[n] -= (a[k]*y[ind])

            k += 1

    return(y)

def filtfilt(b,a,x):

    # filter, flip y, filter again, flip y again
    y = filter(b,a,x)
    y = y.reverse()
    y = filter(b,a,y)
    y = y.reverse()

    return(y)

def iirnotch(w0,bw,Ab=-10*np.log10(1/2)):
    bw = bw * np.pi
    w0 = w0 * np.pi

    Gb = 10**(-Ab/20)
    beta = (np.sqrt(1-Gb**2) / Gb) * np.tan(bw/2)
    gain = 1 / (1+beta)

    num = gain * np.array((1, -2*np.cos(w0), 1))
    den = np.array((1, -2*gain*np.cos(w0), 2*gain-1))

    return(num,den)