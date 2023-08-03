# For this you need matlab and the new requirements.txt
import argparse

import matlab
import matlab.engine
import scipy.io
from oset.ecg.peak_detection.peak_detection_amp_threshold import (
    peak_detection_amp_threshold,
)

import unit_test as testing

mat = scipy.io.loadmat("../../../datasets/sample-data/SampleECG1.mat")["data"][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def peak_detection_amp_threshold_unit_test():
    ml = runMatLab()
    py = runPython()
    return testing.compare_number_arrays(
        py[0], ml[0][0]
    ) and testing.compare_number_arrays(py[1], ml[1][0])


def runMatLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    return eng.peak_detection_amp_threshold(x, f / fs, th, nargout=2)


def runPython():
    return peak_detection_amp_threshold(mat, f / fs, th)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for peak_detection_amp_threshold"""
    )
    args = parser.parse_args()
    print(peak_detection_amp_threshold_unit_test())
