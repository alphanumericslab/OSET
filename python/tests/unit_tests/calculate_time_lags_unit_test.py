import matlab.engine
import scipy.io
import numpy as np
import logging
import argparse
from oset.ecg import phase_calculator
from oset.generic.lp_filter.lp_filter_zero_phase import lp_filter_zero_phase
from oset.ecg.peak_detection.peak_det_simple import peak_det_simple
from oset.ecg.bss.calculate_time_lags import calculate_time_lags
import unit_test as testing

data_path = "../../../datasets/sample-data/SampleECG1.mat"

def calculate_time_lags_unit_test():
    data = scipy.io.loadmat(data_path)["data"][0]
    f = 1
    fs = 1000
    fc = 0.5
    data = data - lp_filter_zero_phase(data, fc / fs)
    peaks, peak_indexes = peak_det_simple(data, f / fs)
    phase, phaseops = phase_calculator(peaks)

    print("Running MATLAB function...")
    t0_mat, t1_mat = run_matlab(peaks, phase)

    print("Running Python function...")
    t0_py, t1_py = run_python(peaks, phase)
    
    print("Comparing outputs...")
    w = testing.compare_number_arrays_with_tolerance(t0_py, t0_mat,tolerance=2)
    x = testing.compare_number_arrays_with_tolerance(t1_py, t1_mat,tolerance=2)
    return w and x

def run_matlab(peaks, phase):
    eng = matlab.engine.start_matlab()
    x = matlab.double(peaks.tolist())
    y = matlab.double(phase.tolist())
    eng.addpath("../../../matlab/tools/ecg/bss")
    t0_mat, t1_mat = eng.calculate_time_lags(x, y, nargout=2)
    eng.quit()
    t0_mat = np.array(t0_mat).flatten()
    t1_mat = np.array(t1_mat).flatten()
    return t0_mat, t1_mat

def run_python(peaks, phase):
    t0_py, t1_py = calculate_time_lags(peaks, phase)
    return t0_py, t1_py

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for calculate_time_lags"""
    )
    args = parser.parse_args()
    print(calculate_time_lags_unit_test())
