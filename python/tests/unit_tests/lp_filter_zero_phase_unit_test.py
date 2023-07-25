# For this you need matlab and the new requirements.txt
import matlab
import matlab.engine
import numpy as np
import scipy.io
import argparse
from oset.generic.lp_filter.lp_filter_zero_phase import lp_filter_zero_phase
import unit_test as testing

mat = scipy.io.loadmat('../../../datasets/sample-data/SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def lp_filter_zero_phase_unit_test():
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for lp_filter_zero_phase"""
    )
    args = parser.parse_args()
    print(lp_filter_zero_phase_unit_test())
