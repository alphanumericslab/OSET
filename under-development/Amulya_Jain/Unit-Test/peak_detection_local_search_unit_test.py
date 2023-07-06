# For this you need matlab and the new requirements.txt
import matlab.engine
import matlab
import scipy.io
from peak_detection_local_search import peak_detection_local_search

mat = scipy.io.loadmat('SampleECG1.mat')['data'][0]
f = 1
fs = 1000


def main():
    ml = runMatLab()
    py = runPython()
    return compare_outputs(py[0], ml[0][0]) and compare_outputs(py[1], ml[1][0])


def runMatLab():
    eng = matlab.engine.start_matlab()
    x = matlab.double(mat.tolist())
    return eng.peak_detection_local_search(x, f / fs, nargout=2)


def runPython():
    return peak_detection_local_search(mat, f / fs)


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
            print(a[i - 2:i + 3])
            print(b[i - 2:i + 3])
            x = False
    return x


if __name__ == "__main__":
    print(main())
