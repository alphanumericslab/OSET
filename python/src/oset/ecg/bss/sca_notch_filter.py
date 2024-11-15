import numpy as np
from scipy.signal import iirnotch, filtfilt
import argparse

def sca_notch_filter(x, fc, Q, fs):
    """
    Multichannel notch filter using Spectral Component Analysis (SCA). Designed
    for removing powerline and its harmonics from multichannel recordings.

    Parameters:
        x (numpy.ndarray): Multichannel input signal (each row represents a channel)
        fc (numpy.ndarray): Array of desired notch frequencies
        Q (numpy.ndarray): Array of desired notch Q factors (in forward path)
        fs (float): Sampling frequency

    Returns:
        tuple: (s, W, A, B, C)
            s (numpy.ndarray): Separated signal after notch filtering
            W (numpy.ndarray): Mixing matrix learned by SCA
            A (numpy.ndarray): Unmixing matrix (inverse of W)
            B (numpy.ndarray): Covariance matrix of the filtered signal
            C (numpy.ndarray): Covariance matrix of the input signal

    Reference:
        Sameni, Reza, Christian Jutten, and Mohammad B. Shamsollahi.
        "A deflation procedure for subspace decomposition." IEEE Transactions on
        Signal Processing 58.4 (2010): 2363-2374.

    Revision History:
        June 2024: Translated to Python from Matlab (sca_notch_filter.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    y = x.copy()

    for k in range(len(fc)):
        wc = fc[k] / (fs / 2)
        bw = wc / Q[k]

        b_notch_filter, a_notch_filter = iirnotch(wc, bw)

        y = filtfilt(b_notch_filter, a_notch_filter, y.T).T

    B = np.cov(y.T)
    C = np.cov(x.T)

    B = (B + B.T) / 2
    C = (C + C.T) / 2

    eigvals, eigvecs = np.linalg.eig(np.linalg.solve(C, B))

    idx = np.argsort(eigvals)[::-1]
    W = eigvecs[:, idx].T
    A = np.linalg.pinv(W)

    s = np.dot(W, x)

    return s, W, A, B, C

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Multichannel notch filter using Spectral Component Analysis (SCA). Designed
    for removing powerline and its harmonics from multichannel recordings.

    Syntax: s, W, A, B, C = sca_notch_filter(x, fc, Q, fs)

    Parameters:
        x (numpy.ndarray): Multichannel input signal (each row represents a channel)
        fc (numpy.ndarray): Array of desired notch frequencies
        Q (numpy.ndarray): Array of desired notch Q factors (in forward path)
        fs (float): Sampling frequency

    Returns:
        tuple: (s, W, A, B, C)
            s (numpy.ndarray): Separated signal after notch filtering
            W (numpy.ndarray): Mixing matrix learned by SCA
            A (numpy.ndarray): Unmixing matrix (inverse of W)
            B (numpy.ndarray): Covariance matrix of the filtered signal
            C (numpy.ndarray): Covariance matrix of the input signal

    Revision History:
        June 2024: Translated to Python from Matlab (sca_notch_filter.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()