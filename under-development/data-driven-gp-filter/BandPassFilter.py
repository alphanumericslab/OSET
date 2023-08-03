# Python implementation for Band Pass filtering a signal

import numpy as np
import scipy.signal as scsig


def BandPassFilter(signal, k_BP, fc):
    alpha = (
        1
        - k_BP * np.cos(2 * np.pi * fc)
        - np.sqrt(
            2 * k_BP * (1 - np.cos(2 * np.pi * fc))
            - k_BP**2 * np.sin(2 * np.pi * fc) ** 2
        )
    ) / (1 - k_BP)
    signal_f = scsig.filtfilt(1 - alpha, [1, -alpha], signal)
    return signal_f
