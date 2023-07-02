# Python function that computes the optimal Tau (i.e. number of samples of the beats in the phase space)
# For every beat in time domain the region between the starting point and the corresponding R-peak is consdiered, and the maximum is selected; same for the region between the R-peak and the ending point; the sum is the optimal Tau;
# Input: s - The ECG
#        beats_Rpeaks - the R-peaks positions (relative to each beat)
# Output: TauFirstHalf - the optimal Tau for the first 'half'
#         TauScondHalf - the optimal Tau for the second 'half'

import numpy as np
def ComputeOptimalTau(s, beats_Rpeaks):
    first_half = []
    second_half = []
    for i in np.arange(len(s)):
        first_half.append(len(s[i][0:beats_Rpeaks[i]]))
        second_half.append(len(s[i][beats_Rpeaks[i]:]))
    TauFirstHalf = max(first_half) 
    TauSecondHalf = max(second_half)
    return TauFirstHalf, TauSecondHalf