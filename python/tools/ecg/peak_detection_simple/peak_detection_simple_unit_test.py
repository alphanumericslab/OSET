# For this you need matlab and the new requirements.txt
import matlab.engine
import matlab
import numpy as np
import scipy.io
from peak_detection_simple import peak_detection_simple

mat = scipy.io.loadmat('../SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def main():
    ml = runMatLab(0)
    py = runPython(0)
    w = compare_outputs(py[0], ml[0][0])
    x = compare_outputs(py[1], ml[1][0])
    del ml, py
    ml = runMatLab(1)
    py = runPython(1)
    y = compare_outputs(py[0], ml[0][0])
    z = compare_outputs(py[1], ml[1][0])
    return w and x and y and z


def runMatLab(z):
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    return eng.peak_detection_simple(x, f / fs, z, nargout=2)


def runPython(z):
    return peak_detection_simple(mat, f / fs, z)


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
        if not (a[i] == b[i]):
            print(i)
            print(a[i - 2:i + 3], "python")
            print(b[i - 2:i + 3], "matlab")
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
