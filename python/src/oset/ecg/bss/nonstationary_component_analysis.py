import argparse
import numpy as np

def nonstationary_component_analysis(x, I, J, mode='COV'):
    """
    nonstationary_component_analysis - Nonstationary Component Analysis (NSCA)
    algorithm to extract temporal nonstationarities from multichannel data

    Parameters:
        x (ndarray): Mixture of signals, where each row represents a different signal (N x T)
        I (array-like): Desired time window indices (in 1:T)
        J (array-like): Desired time window indices (in 1:T)
        mode (str, optional): Calculate covariance ('COV'), correlation ('CORR') or
                              symmetrized covariance ('SYMM-COV') matrices. Default: 'COV'

    Returns:
        y (ndarray): Separated components of the input signals (N x T)
        W (ndarray): Separation filters
        A (ndarray): Inverse of separation filters
        B (ndarray): Covariance matrix computed from time window I
        C (ndarray): Covariance matrix computed from time window J
        lambda (ndarray): Generalized eigenvalues sorted in descending order
        
    Note: The time window indices I and J specify subsets of the input signals
        for performing generalized eigenvalue decomposition. In 'COV' mode the
        time indexes are used to calculate covariance matrices, and in 'CORR' 
        mode correlation matrices are calculated (the channel-wise means are 
        preserved contrasting the 'COV' mode). In 'SYMM-COV' mode the data is
        artificially zero-meaned by mirroring the samples. See reference [3]
        for details of this mode.
        
    Example usage:
        x = ... # Define the mixture of signals
        I = ... # Define the time window indices for covariance matrix B
        J = ... # Define the time window indices for covariance matrix C
        nsca_mode = 'CORR' # Cross correlation mode (preserves channel averages)
        y, W, A, B, C, lambda_ = nonstationary_component_analysis(x, I, J, nsca_mode)

        
    Reference:
        1- Sameni, R., Jutten, C., and Shamsollahi, M. B. (2010). A Deflation
        Procedure for Subspace Decomposition. In IEEE Transactions on Signal
        Processing, (Vol. 58, Issue 4, pp. 2363–2374). doi:
        10.1109/tsp.2009.2037353

        2- Sameni, R., and Gouy-Pailler, C. (2014). An iterative subspace
        denoising algorithm for removing electroencephalogram ocular artifacts.
        In Journal of Neuroscience Methods (Vol. 225, pp. 97–105). doi:
        10.1016/j.jneumeth.2014.01.024
 
        3- Sameni R, Jutten C, Shamsollahi MB. What ICA provides for ECG
        processing: Application to noninvasive fetal ECG extraction. In2006
        IEEE International Symposium on Signal Processing and Information
        Technology 2006 Aug 27 (pp. 656-661). IEEE.

    Revision History:
        June 2024: Translated to Python from Matlab (nonstationary_component_analysis.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    if mode not in ['COV', 'SYMM-COV', 'CORR']:
        raise ValueError('Undefined mode')

    # Compute covariance matrices for the desired time windows
    if mode == 'COV':  # covariance
        B = np.cov(x[:, I], rowvar=True)
        C = np.cov(x[:, J], rowvar=True)
    elif mode == 'SYMM-COV':  # artificially symmetrized covariance. See ref. [3]
        B = np.cov(np.hstack([x[:, I], -x[:, I]]), rowvar=True)
        C = np.cov(np.hstack([x[:, J], -x[:, J]]), rowvar=True)
    elif mode == 'CORR':  # correlation
        B = np.dot(x[:, I], x[:, I].T) / len(I)
        C = np.dot(x[:, J], x[:, J].T) / len(J)

    # Symmetrize the covariance matrices
    C = (C + C.T) / 2
    B = (B + B.T) / 2

    # Perform eigenvalue decomposition using Cholesky decomposition
    eigvals, eigvecs = np.linalg.eigh(np.dot(np.linalg.inv(C), B))

    # Sort eigenvalues in descending order
    sorted_indices = np.argsort(eigvals)[::-1]
    lambda_ = eigvals[sorted_indices]
    V = eigvecs[:, sorted_indices]

    # Extract separation filters
    W = V.T
    A = np.linalg.pinv(W)

    # Apply separation filters to input signals
    y = np.real(np.dot(W, x))

    return y, W, A, B, C, lambda_


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    nonstationary_component_analysis - Nonstationary Component Analysis (NSCA)
    algorithm to extract temporal nonstationarities from multichannel data

    Syntax: y, W, A, B, C, lambda_ = nonstationary_component_analysis(x, I, J, nsca_mode)

    Parameters:
        x (ndarray): Mixture of signals, where each row represents a different signal (N x T)
        I (array-like): Desired time window indices (in 1:T)
        J (array-like): Desired time window indices (in 1:T)
        mode (str, optional): Calculate covariance ('COV'), correlation ('CORR') or
                              symmetrized covariance ('SYMM-COV') matrices. Default: 'COV'

    Returns:
        y (ndarray): Separated components of the input signals (N x T)
        W (ndarray): Separation filters
        A (ndarray): Inverse of separation filters
        B (ndarray): Covariance matrix computed from time window I
        C (ndarray): Covariance matrix computed from time window J
        lambda (ndarray): Generalized eigenvalues sorted in descending order
        
    Revision History:
        June 2024: Translated to Python from Matlab (nonstationary_component_analysis.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()

