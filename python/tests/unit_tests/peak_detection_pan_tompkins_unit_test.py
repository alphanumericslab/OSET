# For this you need matlab and the new requirements.txt
import sys
import matlab.engine
import matlab
import numpy as np
import scipy.io
from peak_detection_pan_tompkins import peak_detection_pan_tompkins
import os

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(module_path)
import Unit_test as testing

mat = scipy.io.loadmat('../../../../datasets/sample-data/SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def main():
    ml = runMatLab()
    py = runPython()
    x = testing.compare_number_arrays(py[0], ml[0][0])
    y = testing.compare_number_arrays(py[1], ml[1][0])
    return x and y


def runMatLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    eng.addpath('../../../../matlab/tools/ecg')
    eng.addpath('../../../../matlab/tools/generic')
    return eng.peak_detection_pan_tompkins(x, np.double(fs), nargout=2)


def runPython():
    return peak_detection_pan_tompkins(mat, fs)


if __name__ == "__main__":
    print(main())
