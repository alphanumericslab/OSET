import numpy as np
import argparse

def phase_calculator(peaks):
    """
    phase_calculator - ECG phase calculation from a given set of R-peaks.

    Usage:
        phase, phasepos = phase_calculator(peaks)

    Inputs:
        peaks: Vector of R-peak pulse train

    Outputs:
        phase: The calculated phases ranging from -pi to pi. The R-peaks are located at phase = 0.
        phasepos: The calculated phases ranging from 0 to 2*pi. The R-peaks are again located at phasepos = 0.

    References:
    - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2008). Multichannel electrocardiogram decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering, 55(8), 1935-1940.
    - Sameni, R., Shamsollahi, M. B., Jutten, C., & Clifford, G. D. (2007). A nonlinear Bayesian filtering framework for ECG denoising. IEEE Transactions on Biomedical Engineering, 54(12), 2172-2185.

    Reza Sameni, 2008-2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    
    Revision History:
        June 2024: Translated to Python from Matlab (phase_calculator.m)

    Muhammad Ubadah Tanveer, 2024
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET

    """
    phasepos = np.zeros(len(peaks))
    I = np.where(peaks)[0]

    for i in range(len(I) - 1):
        m = I[i + 1] - I[i]
        phasepos[I[i] + 1:I[i + 1] + 1] = np.linspace(2 * np.pi / m, 2 * np.pi, m)

    m = I[1] - I[0]
    L = len(phasepos[:I[0] + 1])
    phasepos[:I[0] + 1] = np.linspace(2 * np.pi - (L - 1) * 2 * np.pi / m, 2 * np.pi, L)

    m = I[-1] - I[-2]
    L = len(phasepos[I[-1] + 1:])
    phasepos[I[-1] + 1:] = np.linspace(2 * np.pi / m, L * 2 * np.pi / m, L)

    phasepos = np.mod(phasepos, 2 * np.pi)

    phase = phasepos.copy()
    I = np.where(phasepos > np.pi)[0]
    phase[I] = phasepos[I] - 2 * np.pi

    return phase, phasepos


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    phase_calculator - ECG phase calculation from a given set of R-peaks.

    Syntax: phase, phasepos = phase_calculator(peaks)

    Inputs:
        peaks: Vector of R-peak pulse train

    Outputs:
        phase: The calculated phases ranging from -pi to pi. The R-peaks are located at phase = 0.
        phasepos: The calculated phases ranging from 0 to 2*pi. The R-peaks are again located at phasepos = 0.

    Revision History:
        June 2024: Translated to Python from Matlab (phase_calculator.m)
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()

