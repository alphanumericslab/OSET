# For this you need matlab and the new requirements.txt
import argparse
import matlab
import matlab.engine
import scipy.io
from oset.ecg.peak_detection.peak_detection_matched_filter_robust import peak_detection_matched_filter_robust

import unit_test as testing

mat = scipy.io.loadmat('../../../datasets/sample-data/SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def peak_detection_matched_filter_robust_unit_test():
    import warnings
    warnings.warn("This unit test has not yet been deployed")
    ml = runMatLab(1)
    py = runPython(1)
    a = testing.compare_number_arrays(py[0], ml[0][0])
    b = testing.compare_number_arrays(py[1], ml[1][0])
    ml = runMatLab(2)
    py = runPython(2)
    c = testing.compare_number_arrays(py[0], ml[0][0])
    d = testing.compare_number_arrays(py[1], ml[1][0])
    ml = runMatLab(3)
    py = runPython(3)
    e = testing.compare_number_arrays(py[0], ml[0][0])
    f = testing.compare_number_arrays(py[1], ml[1][0])
    return a and b and c and d and e and f


def runMatLab(itr):
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    y = matlab.double([1, 2, 3, 4])
    eng.addpath('../../../matlab/tools/ecg')
    eng.addpath('../../../matlab/tools/generic')
    return eng.peak_detection_matched_filter_robust(x, f / fs, y, 60, itr, nargout=2)


def runPython(itr):
    return peak_detection_matched_filter_robust(mat, f / fs, [1, 2, 3, 4], 60, itr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for peak_detection_matched_filter_robust"""
    )
    args = parser.parse_args()
    print(peak_detection_matched_filter_robust_unit_test())
