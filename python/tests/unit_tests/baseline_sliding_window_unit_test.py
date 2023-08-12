# For this you need matlab and the new requirements.txt
import argparse

import matlab
import matlab.engine
import scipy.io
from oset.generic.baseline_sliding_window.baseline_sliding_window import (
    baseline_sliding_window,
)

import unit_test as testing

mat = scipy.io.loadmat("../../../datasets/sample-data/SampleECG1.mat")["data"][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def baseline_sliding_window_unit_test():
    ml = runMatLab("md")[0]
    py = runPython("md")
    x = testing.compare_number_arrays(py, ml, round_val=12)
    del ml, py
    ml = runMatLab("mn")[0]
    py = runPython("mn")
    y = testing.compare_number_arrays(py, ml, round_val=12)
    return x and y


def runMatLab(z):
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    temp = matlab.double(50)
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    return eng.baseline_sliding_window(x, temp, z)


def runPython(z):
    return baseline_sliding_window(mat, 50, z)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for baseline_sliding_window"""
    )
    args = parser.parse_args()
    print(baseline_sliding_window_unit_test())
