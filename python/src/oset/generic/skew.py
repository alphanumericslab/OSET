import numpy as np


def skew(data):
    """
    Skewness of a matrix data (over the second dimension).

    Args:
        data (numpy.ndarray): N channels x T samples data matrix.

    Returns:
        tuple: A tuple containing:
            - skw (numpy.ndarray): Skewness vector (N x 1).
            - m (numpy.ndarray): Mean vector (N x 1).
            - sd (numpy.ndarray): Standard deviation vector (N x 1).

    Revision History:
        2024: Translated to Python from Matlab(skew.m)

    Amulya Jain, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    n = data.ndim - 1
    m = np.mean(data, axis=n)  # Compute the mean along the second dimension
    sd = np.std(
        data, axis=n, ddof=1
    )  # Compute the standard deviation along the second dimension
    m3 = np.mean(
        data**3, axis=n
    )  # Compute the mean of data^3 along the second dimension
    # (m3 - 3*m.*sd.^2 - m.^3) ./ sd.^3;
    skw = (m3 - (3 * m * (sd**2)) - (m**3)) / (sd**3)  # Calculate the skewness

    return skw, m, sd
