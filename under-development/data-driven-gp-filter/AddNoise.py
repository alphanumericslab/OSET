# Python implementation for adding Gaussian noise over the signal s for a certain target SNR.
# It returns the noisy signal, the noise itself and the noise variance.

import numpy as np


def AddNoise(target_snr_db, s):
    s = np.asanyarray(s)
    var_n = np.mean(s**2) / (10 ** (target_snr_db / 10))
    n = np.random.normal(0, np.sqrt(var_n), len(s))
    s_n = s + n  # noise up the original signal
    return s_n, n, var_n
