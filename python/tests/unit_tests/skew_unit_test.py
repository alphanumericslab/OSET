# For this you need matlab and the new requirements.txt
import argparse
import matlab
import matlab.engine
import numpy as np
import scipy.io
from oset.generic.skew import skew
import unit_test as testing
import random
import time
from math import isclose


def skew_unit_test():
    random.seed(time.time())
    data = [random.random() for _ in range(1000)]
    data = np.array(data) * 1000
    py = runPython(data)
    ml = runMatLab(data)
    w = isclose(py[0], ml[0], rel_tol=1e-6)
    x = isclose(py[1], ml[1], rel_tol=1e-6)
    y = isclose(py[2], ml[2], rel_tol=1e-6)
    data = [[random.random() for _ in range(1000)] for _ in range(1000)]
    data = np.array(data) * 1000
    py = runPython(data)
    ml = runMatLab(data)
    a = testing.compare_arrays(py[0], ml[0], round_val=6)
    b = testing.compare_arrays(py[1], ml[1], round_val=6)
    c = testing.compare_arrays(py[2], ml[2], round_val=6)
    return w and x and y and a and b and c


def runMatLab(data):
    eng = matlab.engine.start_matlab()
    x = data
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    return eng.skew(x, nargout=3)


def runPython(data):
    return skew(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""This is a unit test for skew""")
    args = parser.parse_args()
    print(skew_unit_test())
