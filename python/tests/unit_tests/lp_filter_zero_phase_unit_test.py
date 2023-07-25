# For this you need matlab and the new requirements.txt
import matlab
import matlab.engine
import numpy as np
import scipy.io
from oset.generic.lp_filter.lp_filter_zero_phase import lp_filter_zero_phase

import unit_test as testing

mat = scipy.io.loadmat('../../../datasets/sample-data/SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def main():
    ml = runMatLab()
    py = runPython()
    x = testing.compare_number_arrays(py, ml[0])
    return x


def runMatLab():
    eng = matlab.engine.start_matlab()
    eng.addpath('../../../matlab/tools/ecg')
    eng.addpath('../../../matlab/tools/generic')
    x = matlab.double(mat.tolist())
    return eng.lp_filter_zero_phase(x, np.double(f / fs), nargout=1)


def runPython():
    return lp_filter_zero_phase(mat, f / fs)


def compare_outputs(a, b):
    x = True
    try:
        if a == b:
            return True
    except:
        print('Iterating through the entire array')
    if not len(a) == len(b):
        raise Exception('lengths of both inputs have to be the same')
    for i in range(len(a)):
        # a[i] = round(a[i], 8)
        # b[i] = round(b[i], 8)
        if not (a[i] == b[i]):
            print(a[i], "Python\n", b[i], "MatLab")
            x = False
    return x


def compare_outputs1(a, b):
    x = True
    try:
        if a == b:
            return True
    except:
        print('Iterating through the entire array')
    if not len(a) == len(b):
        raise Exception('lengths of both inputs have to be the same')
    for i in range(len(a)):
        if not (a[i] == b[i][0]):
            print(i)
            print(a[i - 2:i + 3], "python")
            print(b[i - 2:i + 3], "matlab")
            x = False
    return x


if __name__ == "__main__":
    print(main())
