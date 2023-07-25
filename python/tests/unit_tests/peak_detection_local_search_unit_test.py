# For this you need matlab and the new requirements.txt
import matlab
import matlab.engine
import scipy.io
from oset.ecg.peak_detection.peak_detection_local_search import peak_detection_local_search

import unit_test as testing

mat = scipy.io.loadmat('../../../datasets/sample-data/SampleECG1.mat')['data'][0]
f = 1
fs = 1000


def peak_detection_local_search_unit_test():
    ml = run_matLab()
    py = run_python()
    return testing.compare_number_arrays(py[0], ml[0][0]) and testing.compare_number_arrays(py[1], ml[1][0])


def run_matLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    eng.addpath('../../../matlab/tools/ecg')
    eng.addpath('../../../matlab/tools/generic')
    return eng.peak_detection_local_search(x, f / fs, nargout=2)


def run_python():
    return peak_detection_local_search(mat, f / fs)


if __name__ == "__main__":
    print(peak_detection_local_search_unit_test())
