import argparse

import numpy as np
from scipy.signal import lfilter


def peak_detection_matched_filter(x, fs, h, th, fmax) -> tuple:
    """
    peak_detection_matched_filter - R-peak detector based on a matched filter

    Syntax: [peaks, mn, r] = peak_detection_matched_filter(x, fs, h, th, fmax)

    Args:
        x:      Vector of input data
        fs:     Sampling rate
        h:      Template waveform
        th:     Detection threshold
        fmax:   Maximum expected frequency of the R-peaks

    Returns:
        peaks:  Vector of R-peak impulse train
        mn:     Mean waveform of detected R-peaks
        r:      Filtered output after matched filtering

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_matched_filter.m)

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """
    testing = list()
    n = len(x)
    length = len(h)

    h = np.flip(h)  # Reverse the template waveform

    w = int(np.floor(length / 2))

    r = lfilter(h, 1, np.concatenate((x, np.zeros(w - 1))))  # Matched filtering
    r = r[w - 1 : w + n]  # Trim the filtered output to match the length of x
    r[r < th * max(r)] = 0  # Thresholding: set values below th*max(r) to 0

    peaks = np.zeros(len(x))

    wlen2 = round(fs / fmax)
    i = np.where(r > 0)[0] + 1

    seg = np.zeros((i.shape[0], length))
    for index in range(1, i.shape[0] + 1):
        ind = np.array(
            np.arange(max(1, i[index - 1] - wlen2), min(n, i[index - 1] + wlen2) + 1),
            dtype=np.int64,
        )
        if max(r[ind - 1]) == r[i[index - 1] - 1]:
            peaks[i[index - 1] - 1] = 1
            sg = np.array(
                x[max(1, i[index - 1] - w + 1) - 1 : min(n, i[index - 1] + w)],
                dtype=np.double,
            ).tolist()

            seg[index - 1] = np.concatenate(
                (sg, np.zeros((length - len(sg)) if (length - len(sg)) > 0 else 0))
            )[: seg.shape[1]]

    mn0 = np.mean(seg, 0)
    # Robust weighted averaging
    noise = seg - np.tile(mn0, (seg.shape[0], 1))
    vr = np.var(noise, axis=1)
    sm = np.sum(1 / vr)
    weight = 1 / (vr * sm)
    mn = np.dot(np.transpose(weight), seg)
    mn = np.sqrt(np.sum(np.power(h, 2))) * mn / np.sqrt(np.sum(np.power(mn, 2)))

    return peaks, mn, r


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    peak_detection_matched_filter - R-peak detector based on a matched filter
    
    Syntax: [peaks, mn, r] = peak_detection_matched_filter(x, fs, h, th, fmax)
    
    Args:
        x:      Vector of input data
        fs:     Sampling rate
        h:      Template waveform
        th:     Detection threshold
        fmax:   Maximum expected frequency of the R-peaks
    
    Returns:
        peaks:  Vector of R-peak impulse train
        mn:     Mean waveform of detected R-peaks
        r:      Filtered output after matched filtering

    Revision History:
        July 2023: Translated to Python from Matlab (peak_detection_matched_filter.m)

    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    args = parser.parse_args()
