# Python function that creates tall transformation matrices that project a heartbeat onto the phase domain.
# In order for the model to work, the length of the heartbeats in the phase domain is always greater than the length of the heartbeats in the time domain

# Inputs
#    N - the length of the heartbeat in time domain 
#    M - the length of the heartbeat in the phase domain

# Output
#   mat - the tall M x N transformation matrix that projects a time-domain heartbeat onto the phase-domain heartbeat.

import numpy as np
def ComputeHeartbeatPsi(N, M):
    ratio = N/(M-1)
    output_interval = []
    for i in np.arange(M):
        output_interval.append(i*ratio)
    mat = np.zeros((M,N))    
    for i in np.arange(N):   
        for j in np.arange(M):        
            if (output_interval[j]>=i and output_interval[j]<i+1):
                mat[j,i] = 1
    mat[j,i] = 1
    return mat