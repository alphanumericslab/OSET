import argparse
import numpy as np

def easi_source_separation(x, nsou, lambda_, nlintype):
    """
    easi_source_separation - An adaptive blind source separation algorithm using the EASI (equivariant adaptive separation by independence) algorithm
    
    Parameters:
    x : ndarray
        Input signal as a matrix of size (ncap x T), where ncap is the number of input channels and T is the length of the signal
    nsou : int
        Number of sources
    lambda_ : float
        Adaptation parameter
    nlintype : int
        Nonlinearity type for the adaptation (y2 denotes y**2)
        1: g = np.diag(y2) * y[:, t]
        2: g = y[:, t] * np.diag(0.1 * y2) * np.diag(y2)
        3: g = y[:, t] * np.sqrt(np.diag(y2))
        4: g = y[:, t] * np.log(np.diag(y2))
        5: g = -y[:, t] * (np.diag(y2) < 0.9)
        6: g = y[:, t] / np.log(np.diag(y2))
        7: g = -y[:, t] / np.sqrt(np.diag(y2))
        8: g = -y[:, t] / np.diag(y2)

    Returns:
    y : ndarray
        Separated sources as a matrix of size (nsou x T)
    B : ndarray
        Separation matrix of size (nsou x ncap)
    conv : ndarray
        Convergence of the algorithm at each iteration (1 x T)
    References:
        Beate Laheld and Jean-François Cardoso, "Adaptive source separation without prewhitening," Proc. EUSIPCO'94, 183-186, Edinburgh, Sep. 1994.
   
    Revision History:
        June 2024: Translated to Python from Matlab (easi_source_separation.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    ncap, T = x.shape  # Number of input channels and length of the signal
    idsou = np.eye(nsou)  # Identity matrix for nsou sources

    B = np.random.randn(nsou, ncap)  # Initialization of the separation matrix
    y = np.zeros((nsou, T))  # Initialization of the separated sources
    conv = np.zeros(T)  # Initialization of the convergence vector

    for t in range(T):
        y[:, t] = B @ x[:, t]  # Separation of sources using the separation matrix
        y2 = y[:, t] @ y[:, t]  # Square of the separated sources

        # Nonlinearity adaptation based on the specified nlintype
        if nlintype == 1:
            g = np.diag(y2) @ y[:, t]
        elif nlintype == 2:
            g = y[:, t] @ np.diag(0.1 * y2) @ np.diag(y2)
        elif nlintype == 3:
            g = y[:, t] @ np.sqrt(np.diag(y2))
        elif nlintype == 4:
            g = y[:, t] @ np.log(np.diag(y2))
        elif nlintype == 5:
            g = -y[:, t] @ (np.diag(y2) < 0.9)
        elif nlintype == 6:
            g = y[:, t] / np.log(np.diag(y2))
        elif nlintype == 7:
            g = -y[:, t] / np.sqrt(np.diag(y2))
        elif nlintype == 8:
            g = -y[:, t] / np.diag(y2)

        gy = g @ y[:, t]  # Inner product of the nonlinearity adaptation and the separated source
        G = (y2 - idsou) / (1 + lambda_ * np.trace(y2)) + (gy - gy.T) / (1 + lambda_ * np.abs(g.T @ y[:, t]))  # Update matrix G
        B = B - lambda_ * G @ B  # Update the separation matrix using matrix G
        conv[t] = np.linalg.norm(G)  # Store the convergence value at each iteration

    return y, B, conv

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    easi_source_separation - An adaptive blind source separation algorithm using the EASI (equivariant adaptive separation by independence) algorithm
    
    Parameters:
    x : ndarray
        Input signal as a matrix of size (ncap x T), where ncap is the number of input channels and T is the length of the signal
    nsou : int
        Number of sources
    lambda_ : float
        Adaptation parameter
    nlintype : int
        Nonlinearity type for the adaptation (y2 denotes y**2)
        1: g = np.diag(y2) * y[:, t]
        2: g = y[:, t] * np.diag(0.1 * y2) * np.diag(y2)
        3: g = y[:, t] * np.sqrt(np.diag(y2))
        4: g = y[:, t] * np.log(np.diag(y2))
        5: g = -y[:, t] * (np.diag(y2) < 0.9)
        6: g = y[:, t] / np.log(np.diag(y2))
        7: g = -y[:, t] / np.sqrt(np.diag(y2))
        8: g = -y[:, t] / np.diag(y2)

    Returns:
    y : ndarray
        Separated sources as a matrix of size (nsou x T)
    B : ndarray
        Separation matrix of size (nsou x ncap)
    conv : ndarray
        Convergence of the algorithm at each iteration (1 x T)

    Revision History:
        June 2024: Translated to Python from Matlab (easi_source_separation.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()

