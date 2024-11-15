import numpy as np
import argparse


def chmm_cartesian_product(pi_0_chmm, transition_chmm, chmm_gmm_para):
    """
    Convert a Coupled Hidden Markov Model (CHMM) to an equivalent Hidden Markov Model (EHMM).
    
    Parameters:
    pi_0_chmm : list of np.array
        Initial state probability distributions for each channel in the CHMM.
    transition_chmm : list of np.array
        Transition probability matrices for each channel in the CHMM.
    chmm_gmm_para : list of dict
        Parameters of the Gaussian Mixture Model (GMM) for each state in the CHMM.

    Returns:
    pi_0_ehmm : list of np.array
        Initial state probability distribution for the equivalent HMM.
    coupling_theta_ehmm : int
        Coupling coefficient for the equivalent HMM.
    transition_ehmm : list of np.array
        Transition probability matrix for the equivalent HMM.
    ehmm_gmm_para : list of dict
        Parameters of the GMM for each state in the equivalent HMM.
    index_matrix : np.array
        Mapping from CHMM states to equivalent HMM states.


    Revision History:
        August 2024: Translated to Python from Matlab (chmm_cartesian_product.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    coupling_theta_ehmm = 1

    C = len(pi_0_chmm)
    temp = 1

    for c in range(C):
        temp = np.kron(temp, pi_0_chmm[c])

    pi_0_ehmm = [temp]

    sigma_diag = 1 if chmm_gmm_para[0]['gmm_para'][0]['sigma'][0]['x'].ndim == 1 else 0

    dim_observation = np.zeros(C)
    state_numbers = np.zeros(C, dtype=int)
    num_gmm_component = np.zeros(C, dtype=int)

    for zee in range(C):
        dim_observation[zee] = len(chmm_gmm_para[zee]['gmm_para'][0]['mu'][0]['x'])
        state_numbers[zee] = pi_0_chmm[zee].shape[0]
        num_gmm_component[zee] = len(chmm_gmm_para[zee]['gmm_para'][0]['P'])

    A_cartesian = np.zeros((np.prod(state_numbers), np.prod(state_numbers)))

    for row_num in range(A_cartesian.shape[0]):
        temp = 1
        for c in range(C):
            temp = np.kron(temp, transition_chmm[c][row_num, :])
        A_cartesian[row_num, :] = temp

    transition_ehmm = [A_cartesian]

    index_matrix = np.zeros((C, np.prod(state_numbers)), dtype=int)
    index_matrix_gmm = np.zeros((C, np.prod(num_gmm_component)), dtype=int)

    for c in range(C):
        temp_index = state_numbers.copy()
        temp_index[:c] = []
        temp_raw = np.kron(np.arange(1, state_numbers[c] + 1), np.ones(int(np.prod(temp_index))))
        temp_index = state_numbers.copy()
        temp_index[c:] = []
        index_matrix[c, :] = np.tile(temp_raw, int(np.prod(temp_index)))

        temp_index = num_gmm_component.copy()
        temp_index[:c] = []
        temp_raw = np.kron(np.arange(1, num_gmm_component[c] + 1), np.ones(int(np.prod(temp_index))))
        temp_index = num_gmm_component.copy()
        temp_index[c:] = []
        index_matrix_gmm[c, :] = np.tile(temp_raw, int(np.prod(temp_index)))

    dimension_numbers_index = np.concatenate(([0], np.cumsum(dim_observation).astype(int)))

    channel_num_states = int(np.prod(state_numbers))
    num_gmm_component_hmm = int(np.prod(num_gmm_component))

    P_all_hmm = np.ones((1, channel_num_states, num_gmm_component_hmm))
    mu_all_hmm = np.zeros((sum(dim_observation), channel_num_states, num_gmm_component_hmm))

    if sigma_diag:
        sigma_all_hmm = np.zeros((sum(dim_observation), channel_num_states, num_gmm_component_hmm))
    else:
        sigma_all_hmm = np.zeros((sum(dim_observation), sum(dim_observation), channel_num_states, num_gmm_component_hmm))

    for i in range(channel_num_states):
        this_set = index_matrix[:, i]
        for k in range(num_gmm_component_hmm):
            this_gmm = index_matrix_gmm[:, k]
            for c in range(C):
                mu_all_hmm[dimension_numbers_index[c]:dimension_numbers_index[c + 1], i, k] = (
                    chmm_gmm_para[c]['gmm_para'][this_set[c] - 1]['mu'][this_gmm[c] - 1]['x']
                )
                if sigma_diag:
                    sigma_all_hmm[dimension_numbers_index[c]:dimension_numbers_index[c + 1], i, k] = (
                        chmm_gmm_para[c]['gmm_para'][this_set[c] - 1]['sigma'][this_gmm[c] - 1]['x']
                    )
                else:
                    sigma_all_hmm[dimension_numbers_index[c]:dimension_numbers_index[c + 1],
                                  dimension_numbers_index[c]:dimension_numbers_index[c + 1], i, k] = (
                        chmm_gmm_para[c]['gmm_para'][this_set[c] - 1]['sigma'][this_gmm[c] - 1]['x']
                    )
                P_all_hmm[0, i, k] *= chmm_gmm_para[c]['gmm_para'][this_set[c] - 1]['P'][this_gmm[c] - 1]

    gmm_para = []

    for i in range(channel_num_states):
        gmm_components = {}
        for k in range(num_gmm_component_hmm):
            gmm_components.setdefault('P', []).append(P_all_hmm[0, i, k])
            if sigma_diag:
                gmm_components.setdefault('sigma', []).append({'x': sigma_all_hmm[:, i, k]})
            else:
                gmm_components.setdefault('sigma', []).append({'x': sigma_all_hmm[:, :, i, k]})
            gmm_components.setdefault('mu', []).append({'x': mu_all_hmm[:, i, k]})

        gmm_para.append(gmm_components)

    ehmm_gmm_para = [{'gmm_para': gmm_para}]

    return pi_0_ehmm, coupling_theta_ehmm, transition_ehmm, ehmm_gmm_para, index_matrix


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    Coupled Hidden Markov Model (CHMM) to Equivalent Hidden Markov Model (EHMM) conversion.
    Converts a CHMM with multiple channels to an equivalent single-channel HMM.

    Parameters:
    C (int): Number of channels in the CHMM.
    channel_state_num (list): List of state numbers for each channel.
    channel_obs_dim (list): List of observation dimensions for each channel.
    num_gmm_component (list): List of GMM component numbers for each channel.
    chmm_gmm_para (list): List of CHMM parameters for each channel.
    coupling_theta (ndarray): Coupling matrix between channels.
    sigma_diag (bool): If True, use diagonal covariance matrices. If False, use full covariance matrices.

    Returns:
    tuple: A tuple containing:
        pi_0_ehmm (ndarray): Initial state probabilities for EHMM.
        coupling_theta_ehmm (ndarray): Coupling matrix for EHMM.
        transition_ehmm (ndarray): Transition matrix for EHMM.
        ehmm_gmm_para (list): GMM parameters for EHMM.
        index_matrix (ndarray): Mapping between CHMM and EHMM states.

    Revision History:
        August 2024: Implemented in Python based on MATLAB implementation.
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()