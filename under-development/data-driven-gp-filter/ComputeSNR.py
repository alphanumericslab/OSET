#Python implementation for SNR computation
# Inputs: s - the signal
#         n - the noise
# Inputs: snr - the SNR between s and n

import numpy as np
def ComputeSNR(s, n):
    s = np.asanyarray(s)
    n = np.asanyarray(n)
    snr = 10 * np.log10(np.mean(s**2)/np.mean(n**2))
    return snr

