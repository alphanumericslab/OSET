import numpy as np
from scipy.signal import lfilter
import argparse

def bp_filter_complex_ma(x, fc, bw, order):
    """
    Applies a zero-phase bandpass filter to the input signal using a 
    multi-stage moving average (MA) filter with complex-valued frequency shifting.

    Parameters:
    x (ndarray): Input data vector or matrix (channels x samples).
    fc (float): Normalized center frequency of the bandpass filter (between 0 and 1).
    bw (float): Normalized bandwidth of the bandpass filter (between 0 and 1).
    order (int): Order of the MA filter.

    Returns:
    ndarray: Filtered data vector or matrix (channels x samples).

    Notes:
    - fc and bw are normalized by the sampling frequency.
    - The filter performs forward-reverse filtering successively. It has zero-phase for even MA filter orders 
      and a phase-lag equal to a single stage MA filter for odd MA filter orders.
    
    Revision History:
        July 2024: Translated to Python from Matlab (bp_filter_complex_ma.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    M, N = x.shape  # Number of channels and samples

    # Simple lowpass moving average filter used as a template
    L = round(2 / bw)  # Length of the moving average filter
    h = np.ones(L) / L  # Coefficients of the moving average filter

    # Exponentials used for frequency shifting
    n = np.arange(N)  # Time index
    w = np.exp(-1j * 2 * np.pi * fc * n)  # Complex exponential for forward filtering
    v = np.exp(1j * 2 * np.pi * fc * n)   # Complex exponential for reverse filtering

    # Real part of the filter
    x1 = x * w  # Apply frequency shifting to the input signal

    y1 = x1
    for _ in range(order):
        y1 = np.apply_along_axis(lambda m: np.convolve(m, h, mode='same'), axis=1, arr=y1)
        y1 = np.flip(y1, axis=1)  # Reverse the signal (for reverse filtering)

    if order % 2 == 1:
        y1 = np.flip(y1, axis=1)  # Reverse the signal again for odd filter order

    z1 = y1 * v  # Apply reverse frequency shifting

    # Imaginary part of the filter
    x2 = x * v  # Apply reverse frequency shifting to the input signal

    y2 = x2
    for _ in range(order):
        y2 = np.apply_along_axis(lambda m: np.convolve(m, h, mode='same'), axis=1, arr=y2)
        y2 = np.flip(y2, axis=1)  # Reverse the signal (for reverse filtering)

    if order % 2 == 1:
        y2 = np.flip(y2, axis=1)  # Reverse the signal again for odd filter order

    z2 = y2 * w  # Apply forward frequency shifting

    y = np.real(z1 + z2)  # Combine the real and imaginary parts to obtain the filtered signal

    return y



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Bandpass filter using complex moving average (MA) filtering.
    
    Parameters:
        x (ndarray): Vector or matrix of input data (channels x samples)
        fc (float): Normalized center frequency of the bandpass filter
        bw (float): Normalized bandwidth of the bandpass filter
        order (int): Filter order (number of moving average stages)
        
    Returns:
        ndarray: Vector or matrix of filtered data (channels x samples)
    
    Revision History:
        July 2024: Translated to Python from Matlab (bp_filter_complex_ma.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('x', type=np.ndarray, help='Input data array (channels x samples)')
    parser.add_argument('fc', type=float, help='Normalized center frequency')
    parser.add_argument('bw', type=float, help='Normalized bandwidth')
    parser.add_argument('order', type=int, help='Filter order')
    args = parser.parse_args()