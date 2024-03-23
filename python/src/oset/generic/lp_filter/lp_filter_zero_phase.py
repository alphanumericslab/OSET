import argparse
import numpy as np
from scipy.signal import filtfilt

def lp_filter_zero_phase(x, fc):
    """
    lp_filter_zero_phase - Second-order zero-phase Lowpass filter.

    Syntax: y = lp_filter_zero_phase(x, fc)

    Inputs:
      x: Vector or matrix of input data (channels x samples).
      fc: -3dB cut-off frequency normalized by the sampling frequency.

    Output:
      y: Vector or matrix of filtered data (channels x samples).

      Revision History:
          2006: First release
          2023: Renamed from deprecated version LPFilter()

      Read more:  Mitra, S. (2010). Digital signal processing (4th ed.)
                  New York, NY: McGraw-Hill Professional.

    Revision History:
        2023: Translated to Python from Matlab

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    if fc > 1:
        raise ValueError("fc must be less than 1")
    x = np.array([x], dtype=np.double)
    k = 0.7071  # Cut-off value of 1/sqrt(2) or -6dB amplitude attenuation
    alpha = (
        1
        - k * np.cos(2 * np.pi * fc)
        - np.sqrt(
            (
                2 * k * (1 - np.cos(2 * np.pi * fc))
                - k**2 * np.sin(2 * np.pi * fc) ** 2
            )
        )
    ) / (1 - k)
    y = np.zeros_like(x)
    for i in range(x.shape[0]):
        y[i, :] = filtfilt([1 - alpha], [1, -alpha], x[i, :], padlen=3)
    return y[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    lp_filter_zero_phase - Second-order zero-phase Lowpass filter.

    Syntax: y = lp_filter_zero_phase(x, fc)

    Inputs:
      x: Vector or matrix of input data (channels x samples).
      fc: -3dB cut-off frequency normalized by the sampling frequency.

    Output:
      y: Vector or matrix of filtered data (channels x samples).

      Revision History:
          2006: First release
          2023: Renamed from deprecated version LPFilter()

      Read more:  Mitra, S. (2010). Digital signal processing (4th ed.)
                  New York, NY: McGraw-Hill Professional.

    Revision History:
        2023: Translated to Python from Matlab

    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
