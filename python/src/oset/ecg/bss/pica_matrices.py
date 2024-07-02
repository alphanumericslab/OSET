import argparse
import numpy as np
from oset.ecg import phase_calculator

def pica_matrices(x, peaks):
    """
    Calculates covariance and lagged-covariance matrices required for the pseudo-periodic component analysis algorithm (PiCA).

    Parameters:
        x (np.ndarray): Matrix of input signals (rows represent channels, and columns represent samples)
        peaks (np.ndarray): Vector containing the indices of detected peaks in the input signal

    Returns:
        tuple: A tuple containing:
            - A (np.ndarray): Covariance matrix between the peaks and the lagged time points
            - B (np.ndarray): Covariance matrix of the data

    Reference:
        R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel electrocardiogram decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.

    Revision History:
        June 2024: Translated to Python from Matlab (pica_matrices.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    
    # Calculate the phase and discard the second output (not used)
    phase, _ = phase_calculator(peaks)
    
    # PM time calculation for the peaks sequence
    J = np.where(peaks)[0]
    n1 = np.diff(J)
    prd = round(np.mean(n1))
    wlen = max(n1) - min(n1)
    
    N0 = len(x)
    T1 = np.zeros(len(peaks) - prd - wlen)
    NN = len(T1)
    
    for t in range(NN):
        df = np.abs(phase[t] - phase[max(t + prd - wlen, 0):min(t + prd + wlen, N0)])
        I = np.argmin(df)
        T1[t] = t + prd + I - wlen - 1
    
    T1 = np.maximum(T1, 0)
    T1 = np.minimum(T1, N0 - 1)
    T0 = np.arange(NN)
    
    # Calculate the covariance and lagged-covariance matrices for the peaks sequence
    A = x[:, T0] @ x[:, T1.astype(int)].T
    B = x[:, T0] @ x[:, T0].T
    
    # Symmetrize the covariance matrices
    A = (A + A.T) / 2
    B = (B + B.T) / 2
    
    return A, B


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Calculates covariance and lagged-covariance matrices required for the pseudo-periodic component analysis algorithm (PiCA).

    Parameters:
        x (np.ndarray): Matrix of input signals (rows represent channels, and columns represent samples)
        peaks (np.ndarray): Vector containing the indices of detected peaks in the input signal

    Returns:
        tuple: A tuple containing:
            - A (np.ndarray): Covariance matrix between the peaks and the lagged time points
            - B (np.ndarray): Covariance matrix of the data

    Reference:
        R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel electrocardiogram decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.

    Revision History:
        June 2024: Translated to Python from Matlab (pica_matrices.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()

