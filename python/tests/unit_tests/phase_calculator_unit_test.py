import argparse
import matlab
import matlab.engine
import scipy.io
from oset.ecg import phase_calculator
from oset.generic.lp_filter.lp_filter_zero_phase import lp_filter_zero_phase
from oset.ecg.peak_detection.peak_det_simple import peak_det_simple
import unit_test as testing
import numpy as np

# Load the data
data = scipy.io.loadmat("../../../datasets/sample-data/SampleECG1.mat")["data"][0]

f = 1
fs = 1000
fc = 0.5
t = np.arange(len(data)) / fs
data = data - lp_filter_zero_phase(data, fc / fs)
peaks, peak_indexes = peak_det_simple(data, f / fs)


def phase_calculator_unit_test():
    ml = runMatLab(0)
    py = runPython(0)
    w = testing.compare_number_arrays(py[0], ml[0][0],5,True)
    x = testing.compare_number_arrays(py[1], ml[1][0],5,True)
    del ml, py
    ml = runMatLab(1)
    py = runPython(1)
    y = testing.compare_number_arrays(py[0], ml[0][0],5,True)
    z = testing.compare_number_arrays(py[1], ml[1][0],5,True)
    return w and x and y and z

def runMatLab(z):
    eng = matlab.engine.start_matlab()
    x = matlab.double(peaks.tolist())
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    return eng.phase_calculator(x, nargout=2)

def runPython(z):
    return phase_calculator(peaks)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for phase_calculator"""
    )
    args = parser.parse_args()
    print(phase_calculator_unit_test())
