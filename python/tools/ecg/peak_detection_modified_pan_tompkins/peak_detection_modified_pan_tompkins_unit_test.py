# For this you need matlab and the new requirements.txt
import matlab.engine
import matlab
import numpy as np
import scipy.io
from peak_detection_modified_pan_tompkins import peak_detection_modified_pan_tompkins

mat = scipy.io.loadmat('../SampleECG1.mat')['data'][0]
f = 1
fs = 1000
th = 0.10  # an arbitrary value for testing


def main():
    ml = runMatLab()
    py = runPython()
    x = compare_outputs(py[0], ml[0][0])
    y = compare_outputs1(py[1], ml[1])
    z = compare_outputs1(py[2], ml[2])
    return x and y and z


def runMatLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    return eng.peak_detection_modified_pan_tompkins(x, np.double(fs), nargout=3)


def runPython():
    return peak_detection_modified_pan_tompkins(mat, fs)


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
