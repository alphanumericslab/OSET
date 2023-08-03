import logging
import typing

import numpy as np

"""
    Revision History:
        2023: Created
    
    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
"""


def compare_number_arrays(
    python: typing.Iterable,
    matlab: typing.Iterable,
    round_val: int = -1,
    debug: bool = False,
) -> bool:
    """
    Compares 2 number arrays
    Args:
        python: Python array
        matlab: Matlab array
        round_val (int, optional): Determines whether the inputs should be rounded.
                                   Use -1 for no rounding (default), otherwise specify the number of decimals.
        debug (bool, optional): Enable debug mode. False by default.

    Returns: True if the input are equal, False otherwise.

    """
    if round_val > -1:
        python = np.round(python, round_val)
        matlab = np.round(matlab, round_val)[0]
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    x = True
    try:
        if python == matlab:
            return True
    except:
        logging.debug("Iterating through the entire array")
    if not len(python) == len(matlab):
        print("python length:", len(python))
        print("matlab length:", len(matlab))
        raise Exception("lengths of both inputs have to be the same")
    for i in range(len(python)):
        if not (python[i] == matlab[i]):
            print(i)
            print("python:", python[i])
            print("matlab:", matlab[i])
            x = False
    return x


def compare_arrays(
    python: typing.Iterable,
    matlab: typing.Iterable,
    round_val: int = -1,
    debug: bool = False,
) -> bool:
    """
    Compares an array of numbers ('a') against an array of arrays, where each sub-array contains only one number.
    For example, it compares [1, 2, 3] against [[1], [2], [3]].
    Args:
        python: Python array
        matlab: Matlab array
        round_val (int, optional): Determines whether the inputs should be rounded.
                                   Use -1 for no rounding (default), otherwise specify the number of decimals.
        debug (bool, optional): Enable debug mode. False by default.

    Returns: True if the input are equal, False otherwise.

    """
    if round_val > -1:
        python = np.round(python, round_val)
        for i in range(len(python)):
            matlab[i][0] = round(matlab[i][0], round_val)
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    x = True
    try:
        if python == matlab:
            return True
    except:
        logging.debug("Iterating through the entire array")
    if not len(python) == len(matlab):
        print("python length:", len(python))
        print("matlab length:", len(matlab))
        raise Exception("lengths of both inputs have to be the same")
    for i in range(len(python)):
        if not (python[i] == matlab[i][0]):
            print(i)
            print("python:", python[i])
            print("matlab:", matlab[i][0])
            x = False
    return x


if __name__ == "__main__":
    print("""You can't run this. This is a helper file""")
