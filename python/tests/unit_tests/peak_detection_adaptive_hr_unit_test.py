# For this you need matlab and the new requirements.txt
import argparse

import matlab
import matlab.engine
import scipy.io
from oset.ecg.peak_detection.peak_det_adaptive_hr import (
    peak_det_adaptive_hr,
)

import unit_test as testing

mat = scipy.io.loadmat("../../../datasets/sample-data/SampleECG1.mat")["data"][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def peak_detection_adaptive_hr_unit_test():
    ml = runMatLab(0)
    py = runPython(0)
    x = testing.compare_number_arrays(py, ml[0])
    ml = runMatLab(1)
    py = runPython(1)
    y = testing.compare_number_arrays(py, ml[0])
    return x, y


def runMatLab(z):
    eng = matlab.engine.start_matlab()
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    x = matlab.double(mat.tolist())
    return eng.peak_det_adaptive_hr(
        x, matlab.double(60), matlab.double(fs), z, nargout=2
    )


def runPython(z):
    return peak_det_adaptive_hr(mat, 60, fs, z)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for peak_detection_adaptive_hr"""
    )
    args = parser.parse_args()
    print(peak_detection_adaptive_hr_unit_test())
