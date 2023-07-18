# For this you need matlab and the new requirements.txt
import sys
import matlab.engine
import matlab
import numpy as np
import scipy.io
from peak_detection_modified_pan_tompkins import peak_detection_modified_pan_tompkins

sys.path.append('..')
import Unit_test as testing

mat = scipy.io.loadmat('../SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def peak_detection_modified_pan_tompkins_unit_test():
    ml = run_matLab()
    py = run_python()
    x = testing.compare_number_arrays(py[0], ml[0][0])
    y = testing.compare_arrays(py[1], ml[1])
    z = testing.compare_arrays(py[2], ml[2])
    return x and y and z


def run_matLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    eng.addpath('../../../../matlab/tools/ecg')
    eng.addpath('../../../../matlab/tools/generic')
    return eng.peak_detection_modified_pan_tompkins(x, np.double(fs), nargout=3)


def run_python():
    return peak_detection_modified_pan_tompkins(mat, fs)


if __name__ == "__main__":
    print(peak_detection_modified_pan_tompkins_unit_test())
