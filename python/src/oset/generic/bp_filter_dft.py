import argparse
import numpy as np


import numpy as np

def bp_filter_dft(x, fl, fu):
    """
    Bandpass filter using discrete Fourier transform (DFT) filtering.
    
    Parameters:
        x (ndarray): Vector or matrix of input data (channels x samples).
        fl (float): Normalized lower frequency.
        fu (float): Normalized upper frequency.
        
    Returns:
        ndarray: Vector or matrix of filtered data (channels x samples).
    
    Notes:
    - fl and fu are the lower and upper frequency ranges of the bandpass filter
      normalized by the sampling frequency.
    - The filter does not perform any windowing on the data.
    
    Revision History:
        June 2024: Translated to Python from Matlab (bp_filter_dft.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    
    if fl > 0.5:
        raise ValueError('fl may not exceed 0.5')
    if fu > 0.5:
        raise ValueError('fu may not exceed 0.5')
    if fl > fu:
        raise ValueError('fl may not exceed fu')
    
    N = x.shape[1]
    S = np.fft.fft(x, n=N, axis=1)
    
    k_cut_low = int(np.floor(fl * N))
    S[:, :k_cut_low] = 0
    S[:, (N - k_cut_low + 1):] = 0
    
    k_cut_high = int(np.ceil(fu * N))
    if k_cut_high > 0:
        S[:, (k_cut_high + 1):(N - k_cut_high)] = 0
    else:
        S = 0
    
    y = np.real(np.fft.ifft(S, n=N, axis=1))
    return y

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Bandpass filter using discrete Fourier transform (DFT) filtering
    
    Parameters:
        x (ndarray): Vector or matrix of input data (channels x samples)
        fl (float): Normalized lower frequency.
        fu (float): Normalized upper frequency.
        
    Returns:
        ndarray: Vector or matrix of filtered data (channels x samples)
    
    Revision History:
        June 2024: Translated to Python from Matlab (bp_filter_dft.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
