import numpy as np
from chmm_cartesian_product import chmm_cartesian_product
import argparse


def im_para_eqhmm(pi_0_lsim, lsim_gmm_para, coupling_theta_convex_comb, transition_matrices_convex_comb):
    """
    Generate equivalent HMM parameters from LSIM parameters and CHMM.

    Parameters:
    pi_0_lsim : list of np.array
        Initial state probability distributions for each channel in the LSIM.
    lsim_gmm_para : list of dict
        GMM parameters for each state in the LSIM.
    coupling_theta_convex_comb : np.array
        Coupling coefficients for the convex combination model.
    transition_matrices_convex_comb : list of list of np.array
        Transition matrices for the convex combination model.

    Returns:
    pi_0_ehmm : list of np.array
        Initial state probability distribution for the equivalent HMM.
    coupling_theta_ehmm : int
        Coupling coefficient for the equivalent HMM.
    transition_ehmm : list of np.array
        Transition probability matrix for the equivalent HMM.
    ehmm_gmm_para : list of dict
        GMM parameters for the equivalent HMM.
    index_matrix : np.array
        Mapping from CHMM states to equivalent HMM states.
    pi_0_chmm : list of np.array
        Initial state probability distributions for the CHMM.
    transition_chmm : list of np.array
        Transition matrices for the CHMM.
    chmm_gmm_para : list of dict
        GMM parameters for the CHMM.

    Revision History:
        August 2024: Translated to Python from Matlab (im_para_eqhmm.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET

    """
    C = len(pi_0_lsim)
    channel_num_states = np.zeros(C, dtype=int)

    for c in range(C):
        channel_num_states[c] = len(pi_0_lsim[c])

    # Indexing of Cartesian product
    index_matrix = np.zeros((C, np.prod(channel_num_states)), dtype=int)

    for i in range(C):
        temp_index = channel_num_states.copy()
        temp_index[:i] = []
        temp_raw = np.kron(np.arange(1, channel_num_states[i] + 1), np.ones(int(np.prod(temp_index))))
        temp_index = channel_num_states.copy()
        temp_index[i:] = []
        index_matrix[i, :] = np.tile(temp_raw, int(np.prod(temp_index)))

    weight_state_column = np.cumprod(channel_num_states[::-1][1:])
    weight_state_column = np.concatenate(([1], weight_state_column))[::-1]

    # Generate CHMM transition matrices
    transition_chmm = []
    temp_var = np.zeros(C)

    for zee in range(C):
        transition_chmm_zee = np.zeros((np.prod(channel_num_states), channel_num_states[zee]))

        for i in range(channel_num_states[zee]):
            for j in range(index_matrix.shape[1]):
                column_number_matrix_j = np.dot(weight_state_column, (index_matrix[:, j] - 1)) + 1

                for c in range(C):
                    temp_var[c] = transition_matrices_convex_comb[c][zee][index_matrix[c, j] - 1, i]

                transition_chmm_zee[column_number_matrix_j - 1, i] = np.dot(coupling_theta_convex_comb[:, zee], temp_var)

        transition_chmm.append(transition_chmm_zee)

    pi_0_chmm = pi_0_lsim.copy()
    chmm_gmm_para = lsim_gmm_para.copy()

    # Generating equivalent HMM parameters to perform exact inference
    pi_0_ehmm, coupling_theta_ehmm, transition_ehmm, ehmm_gmm_para = chmm_cartesian_product(
        pi_0_chmm, transition_chmm, chmm_gmm_para
    )

    return pi_0_ehmm, coupling_theta_ehmm, transition_ehmm, ehmm_gmm_para, index_matrix, pi_0_chmm, transition_chmm, chmm_gmm_para


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Interval Markov (IM) to Equivalent Hidden Markov Model (EHMM) conversion.
    Converts an IM model with multiple channels to an equivalent single-channel HMM.

    Parameters:
    pi_0_lsim (list): Initial state probabilities for each channel in the IM model.
    transition_matrices_convex_comb (list): Transition matrices for each channel in the IM model.
    coupling_theta_convex_comb (ndarray): Coupling matrix between channels.
    lsim_gmm_para (list): GMM parameters for each state in the IM model.

    Returns:
    tuple: A tuple containing:
        pi_0_ehmm (ndarray): Initial state probabilities for EHMM.
        coupling_theta_ehmm (ndarray): Coupling matrix for EHMM.
        transition_ehmm (ndarray): Transition matrix for EHMM.
        ehmm_gmm_para (list): GMM parameters for EHMM.
        index_matrix (ndarray): Mapping between IM and EHMM states.
        pi_0_chmm (list): Initial state probabilities for intermediate CHMM.
        transition_chmm (list): Transition matrices for intermediate CHMM.
        chmm_gmm_para (list): GMM parameters for intermediate CHMM.

    Revision History:
        August 2024: Implemented in Python based on MATLAB implementation.
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()