# Python implementation of the GP based filter, USING THE DIAGONAL approximation for the covariances.
# This represents the benchmark (for the DIAGONAL models) - many flavours can be added and potentially improving the results. But this is the implementation that respects the model and offers a benchmark.

# It is the version that is the strict implementation of the theory in the following sense:
    # 1. It has not flavor in terms of using filtered matrices or applying an additional filter on the obtained sample mean in the phase domain.
    # 2. It does respect the thory ad literam in the sense that the Gramian appearing the equations is accounted for. 
# Because it assumes diagonal covariances matrices, the implementiaton was designed to be as optimal as possibe, i.e. working with vectors instead of matrices when possible.

# Inputs: x - the ECG to be filtered
#       : beats_Rpeaks - R-peak positions
#       : TauFirstHalf - Length of the phase-domain beats between the first sample and the R-peak sample
#       : TauSecondHalf - Length of the phase-domain beats between the R-peak sample and the last sample.
#       : fs - sampling frequency
#       : noise_var - the noise variance; used as function parameter. it assumes the estimation procedure to be done outside the function.

# Outputs: sh - the filtred ECG (i.e. the posterior mean corresopnding to the GP)
#          mu_s - the prior mean corresopnding to the GP

import numpy as np
def GaussianProcessFilterInPhaseDiagFast(x, beats_Rpeaks, TauFirstHalf, TauSecondHalf, 
                                  fs, noise_var):
    #------------------------------------------------------------
    from ComputeTransformationPsi import ComputeTransformationPsi
    Ps = ComputeTransformationPsi(x, beats_Rpeaks, TauFirstHalf, TauSecondHalf)
    #------------------------------------------------------------    
    # Compute the phase beat measurements
    xi = []
    for i in np.arange(len(x)):    
        xi.append(np.matmul(Ps[i], x[i]))
    #------------------------------------------------------------
    # Compute the phase beat measurements mean
    mu_xi = np.mean(xi, axis = 0)
    #------------------------------------------------------------
    # Compute the phase beat measurements covariance
    K_xi = np.cov(xi, rowvar = False, bias = False)
    k_xi = np.diag(K_xi)
    #k_si = k_xi - noise_var
    k_si = np.maximum(np.zeros(len(k_xi)),k_xi - noise_var)
    #------------------------------------------------------------
    g = []
    for i in np.arange(len(x)):    
        g.append(np.diag(np.matmul(Ps[i].T, Ps[i])))
    #------------------------------------------------------------
    k_s = []
    for i in np.arange(len(x)):      
        k_s.append(np.matmul(Ps[i].T, k_si) / (g[i]**2))      
    #------------------------------------------------------------ 
    # Compute the time beat prior means (i.e. apply the transoposed matrix transformation on the phase mean)
    mu_s = []
    for i in np.arange(len(x)):
        mu_s.append(np.matmul(Ps[i].T, mu_xi) / g[i])
    #------------------------------------------------------------
    # Compute the time beat filters (i.e. apply the posterior formula using all elements computed above)
    sh = []
    for i in np.arange(len(x)):    
        sh.append(mu_s[i] + k_s[i] / (k_s[i] + noise_var) * (x[i] - mu_s[i]))
    #------------------------------------------------------------
    return sh, mu_s, Ps, g, k_xi, k_si, k_s, mu_xi, xi