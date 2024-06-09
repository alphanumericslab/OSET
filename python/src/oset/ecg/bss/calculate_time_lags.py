import argparse
import numpy as np

def calculate_time_lags(peaks, phase):
    """
    Calculate time lags for the pseudo-periodic component analysis algorithm (PiCA) 
    using the R-peaks and the ECG phase signal.

    Parameters:
        peaks (array-like): Vector containing the indices of detected peaks in the input signal.
        phase (array-like): Vector representing the phase values associated with the peaks (see phase_calculator).

    Returns:
        T0 (numpy.ndarray): Vector containing the time lags for the input peaks.
        T1 (numpy.ndarray): Vector containing the corresponding lagged time points.
        
    Reference:
        R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
        electrocardiogram decomposition using periodic component analysis. IEEE
        Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
        
    Revision History:
        June 2024: Translated to Python from Matlab (calculate_time_lags.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    
    J = np.where(peaks)[0]
    n1 = np.diff(J)
    prd = round(np.mean(n1))
    wlen = max(n1) - min(n1)

    T1 = np.zeros(len(peaks) - prd - wlen)
    NN = len(T1)
    
    for t in range(NN):
        df = np.abs(phase[t] - phase[t + prd - wlen: t + prd + wlen + 1])
        I = np.argmin(df)
        T1[t] = t + prd + I - wlen - 1

    T0 = np.arange(1, NN + 1)
    
    return T0, T1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Calculate time lags for the pseudo-periodic component analysis algorithm (PiCA) 
    using the R-peaks and the ECG phase signal.

    Parameters:
        peaks (array-like): Vector containing the indices of detected peaks in the input signal.
        phase (array-like): Vector representing the phase values associated with the peaks (see phase_calculator).

    Returns:
        T0 (numpy.ndarray): Vector containing the time lags for the input peaks.
        T1 (numpy.ndarray): Vector containing the corresponding lagged time points.
        
    Revision History:
        June 2024: Translated to Python from Matlab (calculate_time_lags.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()

