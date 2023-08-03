# BASELINE WANDER EXTRACTION
# function for baseline wander extraction. Based on a moving average (MA) and median (MD) filter.
# the original signal is filtered with a MA filter, then with a MD filter, obtaining the baselinewonder.
# the filtered signal is the difference between the original signal and the double filtered signal.
# Inputs:
# The signal, the MA windows length (win_l_MA), the MD window's length (win_l_MA) and the signal's sampling frequency.
# Returns:
# The baselinewander (bw) and the signal with bw removed (signal_bwr)
#
import numpy as np
from BandPassFilter import BandPassFilter


def BaselineWanderRemovalBP(signal, k_BP, fc1, fc2):
    signal_bwr = signal - BandPassFilter(signal, k_BP, fc1)
    signal_bwr = BandPassFilter(signal_bwr, k_BP, fc2)
    return signal_bwr
