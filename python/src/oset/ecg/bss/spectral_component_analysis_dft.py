import numpy as np
from scipy.linalg import eigh
import argparse

def spectral_component_analysis_dft(x, fl, fu):
    """
    Spectral component analysis (SCA) using a Fourier domain approach.
    Extracts and ranks linear mixtures of multichannel data with maximal energy
    in a given frequency band.

    Parameters:
    x (ndarray): Input data array (channels x samples).
    fl (float): Lower cutoff frequency normalized by sampling frequency.
    fu (float): Upper cutoff frequency normalized by sampling frequency.

    Returns:
    tuple: A tuple containing:
        y (ndarray): Extracted spectral components ranked by their energy in the frequency band of interest.
        W (ndarray): Extraction (separation) matrix.
        A (ndarray): Mixing matrix.

    Reference:
        Sameni, R., Jutten, C., & Shamsollahi, M. B. (2009). A deflation procedure for subspace decomposition. 
        IEEE Transactions on Signal Processing, 58(4), 2363-2374.
    
    Revision History:
        July 2024: Translated to Python from Matlab (spectral_component_analysis_dft.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    L2 = x.shape[1]

    freq_indexes = np.arange(max(round(L2 * fl), 1), min(round(L2 * fu), L2) + 1)

    X = np.fft.fft(x, L2, axis=1)
    Bf = X[:, freq_indexes] @ X[:, freq_indexes].conj().T
    B = 2 * np.real(Bf)
    B = (B + B.T) / 2

    C = x @ x.T
    C = (C + C.T) / 2

    D, V = eigh(B, C)
    II = np.argsort(D)[::-1]

    W = V[:, II].T
    A = np.linalg.pinv(W)

    y = np.real(W @ x)

    return y, W, A


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Spectral component analysis (SCA) using a Fourier domain approach.
    Extracts and ranks linear mixtures of multichannel data with maximal energy
    in a given frequency band.

    Parameters:
    x (ndarray): Input data array (channels x samples).
    fl (float): Lower cutoff frequency normalized by sampling frequency.
    fu (float): Upper cutoff frequency normalized by sampling frequency.

    Returns:
    tuple: A tuple containing:
        y (ndarray): Extracted spectral components ranked by their energy in the frequency band of interest.
        W (ndarray): Extraction (separation) matrix.
        A (ndarray): Mixing matrix.

    Revision History:
        June 2024: Translated to Python from Matlab (spectral_component_analysis_dft.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()