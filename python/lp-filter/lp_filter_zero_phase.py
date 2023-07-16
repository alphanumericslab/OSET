import numpy as np
from scipy.signal import filtfilt

__author__ = 'Reza Sameni'
__source__ = 'The Open-Source Electrophysiological Toolbox (OSET) - https://github.com/alphanumericslab/OSET'
__reference__ = 'Mitra, S. (2010). Digital signal processing (4th ed.) New York, NY: McGraw-Hill Professional.'

def lp_filter_zero_phase(x, fc):
    """
    lp_filter_zero_phase - Second-order zero-phase Lowpass filter.

    Args:
        x: Vector or matrix of input data (channels x samples).
        fc: -3dB cut-off frequency normalized by the sampling frequency.

    Returns:
        y: Vector or matrix of filtered data (channels x samples).
    """
    k = 0.7071  # Cut-off value of 1/sqrt(2) or -6dB amplitude attenuation
    alpha = (1 - k * np.cos(2 * np.pi * fc) - np.sqrt(2 * k * (1 - np.cos(2 * np.pi * fc)) - k ** 2 * np.sin(2 * np.pi * fc) ** 2)) / (1 - k)
    y = np.zeros_like(x)

    for i in range(x.shape[0]):
        y[i, :] = filtfilt([1 - alpha], [1, -alpha], x[i, :])

    return y

