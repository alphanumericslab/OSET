import numpy as np


def tanh_saturation(x, ksigma):
    """
    Saturates outlier samples using a tanh shape function.

    Args:
        x (ndarray): Input data, can be a vector or a matrix (channels x time).
        ksigma (float): Scaling factor for the saturation level.

    Returns:
        ndarray: Saturated data with outliers replaced by the saturation level.

    References:
        Reza Sameni, 2020-2023
        The Open-Source Electrophysiological Toolbox
        https://github.com/alphanumericslab/OSET
    """
    # Compute the scaling factor based on the standard deviation of each channel
    alpha = ksigma * np.std(x, axis=0, keepdims=True)

    # Convert to a column vector
    alpha = alpha.flatten()

    # Scale the input data and apply the tanh function to saturate outliers
    y = (np.diag(alpha) * np.tanh(np.diag(1. / alpha) * x))[0]
    return y
