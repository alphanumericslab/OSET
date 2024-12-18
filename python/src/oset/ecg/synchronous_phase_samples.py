import argparse
import numpy as np

def synchronous_phase_samples(peaks, phase=None):
    """
    synchronous_phase_samples - Calculation of synchronous time instants from beat to
    beat given a set of R-peaks, required for the periodic component analysis
    (PiCA) algorithm.

    Syntax:
      T0, T1 = synchronous_phase_samples(peaks, phase)

    Parameters:
      peaks: Vector of R-peak pulse train.
      phase (optional): The ECG phase obtained from phase_calculator

    Returns:
      T0: First (or reference) time instant vector.
      T1: Second time vector having synchronous phases with T0.

    Reference:
    - R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel electrocardiogram
      decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering,
      55(8):1935-1940, Aug. 2008.

    Revision History:
      June 2024: Translated to Python from Matlab (synchronous_phase_samples.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    I_peaks = np.where(peaks)[0]
    D = np.diff(I_peaks)
    
    if phase is not None:  # Use pre-calculated ECG phase
        L = len(peaks)
        prd = round(np.mean(D))
        wlen = max(D) - min(D)

        T1 = np.zeros(L - prd - wlen)
        NN = len(T1)
        for t in range(NN):
            df = np.abs(phase[t] - phase[max(t + prd - wlen, 0): min(t + prd + wlen, L)])
            I = np.argmin(df)
            T1[t] = t + prd + I - wlen - 1
        
        T1 = np.maximum(T1, 1)
        T1 = np.minimum(T1, NN)
        T0 = np.arange(1, NN + 1)
    
    else:
        if len(I_peaks) < 3:
            T0 = []
            T1 = []
        else:
            start = I_peaks[0]
            stop = I_peaks[-2]

            T1 = np.zeros(stop - start + 1)
            k = 0
            for t in range(start, stop + 1):
                T1[t - start] = I_peaks[k + 1] + round((t - I_peaks[k]) * D[k + 1] / D[k])
                if t >= I_peaks[k + 1]:
                    k += 1
            T0 = np.arange(start, stop + 1)
    
    return T0, T1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    synchronous_phase_samples - Calculation of synchronous time instants from beat to
    beat given a set of R-peaks, required for the periodic component analysis
    (PiCA) algorithm.

    Syntax:
      T0, T1 = synchronous_phase_samples(peaks, phase)

    Parameters:
      peaks: Vector of R-peak pulse train.
      phase (optional): The ECG phase obtained from phase_calculator

    Returns:
      T0: First (or reference) time instant vector.
      T1: Second time vector having synchronous phases with T0.

    Revision History:
      June 2024: Translated to Python from Matlab (synchronous_phase_samples.m)

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
