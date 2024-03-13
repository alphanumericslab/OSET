import numpy as np
from scipy.linalg import toeplitz

def notch_filter_cls(x, ff, gamma):
    """
    A constrained least squares notch filter.

    Apply a constrained least squares notch filter to an input signal.
    This filter is designed to cancel mains frequency noise from a signal.

    Parameters:
    x (numpy.ndarray): Input signal (channels x time).
    ff (float): Mains frequency divided by the sampling frequency.
    gamma (float): Regularization factor.

    Returns:
    numpy.ndarray: The denoised signal (after mains cancellation).

    Method:
    The powerline signal s_k is known to satisfy:
        I) s_{k+1} + s_{k-1} = 2*cos(2*pi*f_mains/fs) * s_k.
    Signals corrupted by the powerline noise can be modeled as:
        II) x_k = s_k + n_k
    (I) and (II) can be combined to solve the constrained least squares
    problem per channel: S_opt = argmin(|H*S| + gamma*|X - S|) where H is
    an oscillator's equation in Toeplitz form (from I), X is the vector form
    of the samples x_k and S is the vector form of s_k. The norm represents 
    the Frobenius norm. The solution is known to be \hat{S} = (I + gamma*H^T*H)^{-1} * X

    NOTE: Due to the large matrix inversion involved, this method is only
    computationally efficient for relatively short-lengthed signals.
    Consider using other notch filter methods for long records.
    """
    N = x.shape[1]
    col = np.concatenate(([1], np.zeros(N - 3)))
    row = np.concatenate(([1], [-2 * np.cos(2 * np.pi * ff)], [1], np.zeros(N - 3)))
    H = toeplitz(col, row)
    A = np.eye(N) + gamma * np.matmul(H.T, H)
    y = x - np.linalg.solve(A, x.T).T

    return y
