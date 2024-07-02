import matlab.engine
import scipy.io
import numpy as np
import argparse
from oset.generic.lp_filter.lp_filter_zero_phase import lp_filter_zero_phase
from oset.ecg.peak_detection.peak_det_simple import peak_det_simple
import unit_test as testing

data_path = "../../../datasets/sample-data/SampleECG1.mat"

def synchronous_phase_samples_unit_test():
    data = scipy.io.loadmat(data_path)["data"][0]
    f = 1
    fs = 1000
    fc = 0.5
    data = data - lp_filter_zero_phase(data, fc / fs)
    peaks, peak_indexes = peak_det_simple(data, f / fs)
    from oset.ecg.phase_calculator import phase_calculator
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
    eng.addpath("../../../matlab/tools/ecg")
    t0_mat, t1_mat = eng.synchronous_phase_samples(x, y, nargout=2)
    eng.quit()
    t0_mat = np.array(t0_mat).flatten()
    t1_mat = np.array(t1_mat).flatten()
    return t0_mat, t1_mat

def run_python(peaks, phase):
    from oset.ecg import synchronous_phase_samples
    t0_py, t1_py = synchronous_phase_samples(peaks, phase)
    return t0_py, t1_py

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This is a unit test for synchronous_phase_samples"""
    )
    args = parser.parse_args()
    print(synchronous_phase_samples_unit_test())
