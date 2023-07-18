import numpy as np
from scipy.signal import lfilter
from tanh_saturation import tanh_saturation
from lp_filter_zero_phase import lp_filter_zero_phase


def peak_detection_modified_pan_tompkins(data, fs, *args):
    """
    peak_detection_modified_pan_tompkins - R-peak detector based on modified
      Pan-Tompkins method. The filters and post-detection R-peak selection
      logic differ from the original algorithm

      peaks, peak_indexes, width = peak_detection_modified_pan_tompkins(data, fs, wlen, fp1, fp2, th, ksigma, flag)

      Inputs:
          data: Vector of input data
          fs: Sampling rate
          wlen: Optional. Moving average window length (default = 150ms)
          fp1: Optional. Lower cut-off frequency (default = 10Hz)
          fp2: Optional. Upper cut-off frequency (default = 33.3Hz)
          th: Optional. Detection threshold (default = 0.2)
          ksigma: Optional. Saturates peaks at ksigma x STD of the signal after energy envelope calculation (default = 12)
          flag: Optional. Search for positive (flag=1) or negative (flag=0) peaks.
                By default, the maximum absolute value of the signal determines the peak sign.

      Outputs:
          peaks: Vector of R-peak impulse train
          peak_indexes: Vector of R-peak indexes
          width: Rise to fall width of the signal's peak bump (in samples)

      Reference:
          Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
          Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532

    Revision History:
        2023: Translated to Python from Matlab

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """

    nargin = len(args) + 2

    if nargin > 2 and args[0] is not None:
        wlen = args[0]
    else:
        wlen = 0.150  # moving average window length 150ms

    if nargin > 3 and args[1] is not None:
        fp1 = args[1]
    else:
        fp1 = 10  # First zero of the HP filter is placed at f = 10Hz

    if nargin > 4 and args[2] is not None:
        fp2 = args[2]
    else:
        fp2 = 33.3  # First zero of the LP filter is placed at f = 33.3Hz

    if nargin > 5 and args[3] is not None:
        th = args[3]
    else:
        th = 0.2

    if nargin > 6 and args[4] is not None:
        flag = args[4]
    else:
        flag = np.abs(np.max(data)) > np.abs(np.min(data))
    if nargin > 7 and args[5] is not None:
        ksigma = args[5]
    else:
        ksigma = 12

    N = len(data)
    data = np.array(data).flatten()

    L1 = round(fs / fp2)
    L2 = round(fs / fp1)

    x0 = data - lp_filter_zero_phase(data, 0.05 / fs)

    x = lfilter(np.concatenate(([1], np.zeros(L1 - 1), [-1])), [L1, -L1], x0)
    x = lfilter(np.concatenate(([1], np.zeros(L1 - 1), [-1])), [L1, -L1], x)
    x = np.concatenate((x[L1 - 1:], np.zeros(L1 - 1) + x[-1]))  # Lag compensation

    # HP filter
    y = lfilter(np.concatenate(([L2 - 1, -L2], np.zeros(L2 - 2), [1])), [L2, -L2], x)

    # Differentiation
    z = np.diff(np.concatenate(([y[0]], y)))

    # Squaring
    w = np.square(z)

    # Moving average
    L3 = round(fs * wlen)
    v = lfilter(np.concatenate(([1], np.zeros(L3 - 1), [-1])), [L3, -L3], w)
    v = np.concatenate([v[round(L3 / 2) - 1:], np.zeros(round(L3 / 2) - 1) + v[-1]])  # Group-delay lag compensation

    v_sat = tanh_saturation(v, ksigma)

    p = v_sat > (th * np.max(v_sat))

    # Edge detection
    rising = np.where(np.diff(np.concatenate(([0], p))) == 1)[0] + 1
    falling = np.where(np.diff(np.concatenate((p, [0]))) == -1)[0] + 1

    if len(rising) == len(falling) - 1:
        rising = np.concatenate(([1], rising))
    elif len(rising) == len(falling) + 1:
        falling = np.concatenate((falling, [N]))

    peak_indexes = np.zeros(len(rising))
    width = np.zeros(len(rising))

    if flag:
        for i in range(len(rising)):
            mx = np.argmax(data[rising[i] - 1:falling[i]])
            peak_indexes[i] = mx - 1 + rising[i]
            width[i] = falling[i] - rising[i]
    else:
        for i in range(len(rising)):
            mn = np.argmin(data[rising[i] - 1:falling[i]])
            peak_indexes[i] = mn - 1 + rising[i]
            width[i] = falling[i] - rising[i]

    peaks = np.zeros(N)
    peaks[peak_indexes.astype(int)] = 1
    peak_indexes += 1
    return peaks, peak_indexes, width
