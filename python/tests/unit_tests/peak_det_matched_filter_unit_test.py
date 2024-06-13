# For this you need matlab and the new requirements.txt
import argparse

import matlab
import matlab.engine
import scipy.io
from oset.ecg.peak_detection.peak_det_matched_filter import (
    peak_det_matched_filter,
)

import unit_test as testing

mat = scipy.io.loadmat("../../../datasets/sample-data/SampleECG1.mat")["data"][0]
f = 1
fs = 1000
th = 0.90  # an arbitrary value for testing


def peak_det_matched_filter_unit_test():
    import warnings

    ml = runMatLab()
    py = runPython()
    x = testing.compare_number_arrays(py[0], ml[0][0])
    y = testing.compare_number_arrays(py[1], ml[1][0], round_val=12)
    z = testing.compare_number_arrays(py[2], ml[2][0])
    return x and y and z


def runMatLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    y = matlab.double([1, 2, 3, 4])
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    return eng.peak_det_matched_filter(x, f / fs, y, th, 60, nargout=3)


def runPython():
    return peak_det_matched_filter(mat, f / fs, [1, 2, 3, 4], th, 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for peak_det_matched_filter"""
    )
    args = parser.parse_args()
    print(peak_det_matched_filter_unit_test())
