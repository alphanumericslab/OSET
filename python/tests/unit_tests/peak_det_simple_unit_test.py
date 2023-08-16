# For this you need matlab and the new requirements.txt
import argparse

import matlab
import matlab.engine
import scipy.io
from oset.ecg.peak_detection.peak_det_simple import peak_det_simple

import unit_test as testing

mat = scipy.io.loadmat("../../../datasets/sample-data/SampleECG1.mat")["data"][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def peak_det_simple_unit_test():
    ml = runMatLab(0)
    py = runPython(0)
    w = testing.compare_number_arrays(py[0], ml[0][0])
    x = testing.compare_number_arrays(py[1], ml[1][0])
    del ml, py
    ml = runMatLab(1)
    py = runPython(1)
    y = testing.compare_number_arrays(py[0], ml[0][0])
    z = testing.compare_number_arrays(py[1], ml[1][0])
    return w and x and y and z


def runMatLab(z):
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    return eng.peak_det_simple(x, f / fs, z, nargout=2)


def runPython(z):
    return peak_det_simple(mat, f / fs, z)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for peak_det_simple"""
    )
    args = parser.parse_args()
    print(peak_det_simple_unit_test())
