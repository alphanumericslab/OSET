# Python function that computes the 
#    input, output and delta SNR 

# Inputs:
#     the (input) signal x
#     the noise (i.e. the noise contained in x) n 
#     the (output, i.e. the filtered signal) s

# Ouputs:
#     input_SNR_db - the SNR computed for the (input) signal x and the noise n
#     output_SNR_db - the SNR computed for the (input) signal x and the "estimated" noise, i.e. x-s
#     delta_SNR_db - difference between output_SNR_db and input_SNR_db
    
from ComputeSNR import ComputeSNR
def ComputeInputOutputDeltaSNRs(x, n, s):
    input_SNR_db = ComputeSNR(x, n)
    output_SNR_db = ComputeSNR(x, [s_elem - x_elem for (s_elem, x_elem) in zip(s, x)])
    delta_SNR_db = output_SNR_db - input_SNR_db
    return input_SNR_db, output_SNR_db, delta_SNR_db

