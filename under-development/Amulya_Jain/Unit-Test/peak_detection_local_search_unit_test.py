# for this you need matlab and the new requirements.txt
import matlab.engine
import matlab
import scipy.io
import numpy as np
import peak_detection_local_search

mat = scipy.io.loadmat('SampleECG1.mat')['data'][0]
f = 1
fs = 1000


def main():
    mlold = runMatLabold()
    mlnew = runMatLaboldnew()
    py = runPython()
    print(compare_outputs(py[0], mlnew[0][0]))
    #print(compare_outputs(py[0], mlold[0][0]))
    #print(compare_outputs(mlold[0][0], mlnew[0][0]))


def runMatLaboldnew():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    return eng.peak_detection_local_search(x, f / fs, nargout=2)


def runMatLabold():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    return eng.PeakDetection(x, f / fs, nargout=2)


def runPython():
    return peak_detection_local_search.peak_detection(mat, f / fs)


def compare_outputs(a, b):
    x = True
    try:
        if a == b:
            return True
    except:
        print('iterating through the entire array')
    if not len(a) == len(b):
        raise Exception('lengths of both inputs have to be the same')
    for i in range(len(a)):
        if not (a[i] == b[i]):
            print(i)
            print(a[i - 2:i + 3])
            print(b[i - 2:i + 3])
            x = False
    return x


if __name__ == "__main__":
    main()
