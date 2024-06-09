import argparse
import numpy as np

def mica_projectors(A):
    """
    mica_projectors - Multidimensional Independent Component Analysis (MICA) projection matrices.

    Description:
        Given an nxn matrix A, which has been estimated by, for example,
        independent component analysis of mixtures following the model x = A*s,
        this function provides the projectors that can be used to decompose x
        into its linear decomposition as follows: x = \sigma_p x_p. Using the
        projectors we can find: x_p = P_tilde[p] @ x. See the reference below
        for details and notations.

    Parameters:
        A : ndarray
            An nxn matrix representing the estimated mixing matrix in the ICA model.

    Returns:
        P_tilde : list of ndarray
            A list of n projection matrices (each n x n) for each independent component.
        P : list of ndarray
            A list of n projection matrices (each n x n) used to create P_tilde.

    Reference:
        Cardoso, J-F. "Multidimensional independent component analysis." In
        Proceedings of the 1998 IEEE International Conference on Acoustics,
        Speech and Signal Processing, ICASSP'98 (Cat. No. 98CH36181), vol. 4,
        pp. 1941-1944. IEEE, 1998. doi: 10.1109/ICASSP.1998.681443

    Revision History:
        June 2024: Translated to Python from Matlab (mica_projectors.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    n = A.shape[1]  # Get the number of columns in A (assumes square matrix).

    P = [None] * n  # Initialize a list for P.
    S = np.zeros((n, n))  # Initialize S as an n x n matrix of zeros.

    for k in range(n):
        # Calculate the projection matrix P[k] for each independent component.
        P[k] = np.outer(A[:, k], A[:, k]) / (A[:, k].T @ A[:, k])
        S += P[k]  # Accumulate P[k] into S.

    P_tilde = [None] * n  # Initialize a list for P_tilde.
    S_inv = np.linalg.pinv(S)  # Calculate the pseudo-inverse of S.

    for k in range(n):
        # Calculate P_tilde by multiplying each P[k] with S_inv.
        P_tilde[k] = P[k] @ S_inv

    return P_tilde, P

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    mica_projectors - Multidimensional Independent Component Analysis (MICA) projection matrices.

    Syntax: P_tilde, P = mica_projectors(A)
    
    Parameters:
        A : ndarray
            An nxn matrix representing the estimated mixing matrix in the ICA model.

    Returns:
        P_tilde : list of ndarray
            A list of n projection matrices (each n x n) for each independent component.
        P : list of ndarray
            A list of n projection matrices (each n x n) used to create P_tilde.

    Revision History:
        June 2024: Translated to Python from Matlab (mica_projectors.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()

