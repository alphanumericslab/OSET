import matlab.engine
import scipy.io
import numpy as np
import logging
import typing
import argparse
from oset.ecg import phase_calculator
from oset.generic.lp_filter.lp_filter_zero_phase import lp_filter_zero_phase
from oset.ecg.peak_detection.peak_det_simple import peak_det_simple


def compare_arrays_with_tolerance(
    python: typing.Iterable,
    matlab: typing.Iterable,
    decimal_places: int = -1,
    debug: bool = False,
) -> bool:
    """
    Compares 2 number arrays with specified tolerance.

    Args:
        python: Python array.
        matlab: Matlab array.
        decimal_places (int, optional): Determines whether the inputs should be compared up to certain decimals.
                                        Use -1 for exact comparison (default), otherwise specify the number of decimals.
        debug (bool, optional): Enable debug mode. False by default.

    Returns: True if the inputs are equal up to the specified decimal places, False otherwise.
    """
    if debug:
        logging.basicConfig(level=logging.DEBUG, force=True)

    if decimal_places > -1:
        atol = 10 ** (-decimal_places)
    else:
        atol = 0

    if len(python) != len(matlab):
        if debug:
            logging.debug(f"Python length: {len(python)}, Matlab length: {len(matlab)}")
        raise ValueError("Lengths of both inputs must be the same")

    comparison_result = True
    discrepancies = []

    for i in range(len(python)):
        if not np.isclose(python[i], matlab[i], atol=atol):
            if (np.isclose(python[i], np.pi, atol=atol) and np.isclose(matlab[i], -np.pi, atol=atol)) or \
               (np.isclose(python[i], -np.pi, atol=atol) and np.isclose(matlab[i], np.pi, atol=atol)):
                logging.debug(f"Index {i}: Special case for pi and -pi handled")
                continue
            discrepancies.append((i, python[i], matlab[i]))
            logging.debug(f"Index {i}: Python {python[i]}, Matlab {matlab[i]}")
            comparison_result = False

    if discrepancies:
        logging.debug(f"Total discrepancies: {len(discrepancies)}")
        for idx, py_val, mat_val in discrepancies:
            logging.debug(f"Index {idx}: Python {py_val}, Matlab {mat_val}")

    return comparison_result


def phase_calculator_unit_test(data_path: str) -> bool:
    data = scipy.io.loadmat(data_path)["data"][0]
    f = 1
    fs = 1000
    fc = 0.5
    t = np.arange(len(data)) / fs
    data = data - lp_filter_zero_phase(data, fc / fs)
    peaks, peak_indexes = peak_det_simple(data, f / fs)

    ml = run_matlab(peaks)
    py = run_python(peaks)
    w = compare_arrays_with_tolerance(py[0], ml[0][0], decimal_places=5, debug=True)
    x = compare_arrays_with_tolerance(py[1], ml[1][0], decimal_places=5, debug=True)
    return w and x


def run_matlab(peaks) -> typing.Tuple[np.ndarray, np.ndarray]:
    eng = matlab.engine.start_matlab()
    x = matlab.double(peaks.tolist())
    eng.addpath("../../../matlab/tools/ecg")
    eng.addpath("../../../matlab/tools/generic")
    result = eng.phase_calculator(x, nargout=2)
    eng.quit()
    return result


def run_python(peaks) -> typing.Tuple[np.ndarray, np.ndarray]:
    return phase_calculator(peaks)


if __name__ == "__main__":
    data_path = "../../../datasets/sample-data/SampleECG1.mat"
    result = phase_calculator_unit_test(data_path)
    print("Test result:", result)