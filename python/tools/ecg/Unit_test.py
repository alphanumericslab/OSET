import logging

import numpy as np

"""
    Revision History:
        2023: Created
    
    Amulya Jain, 2023
    The Open-Source Electrophysiological Toolbox
    https://github.com/alphanumericslab/OSET
"""


def compare_number_arrays(a, b, round_val: int = -1, debug: bool = False) -> bool:
    """
    Compares 2 number arrays
    Args:
        a: Python array
        b: Matlab array
        round_val (int, optional): Determines whether the inputs should be rounded.
                                   Use -1 for no rounding (default), otherwise specify the number of decimals.
        debug (bool, optional): Enable debug mode. False by default.

    Returns: True if the input are equal, False otherwise.

    """
    if round_val > -1:
        a = np.round(a, round_val)
        b = np.round(b, round_val)
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    x = True
    try:
        if a == b:
            return True
    except:
        logging.debug('Iterating through the entire array')
    if not len(a) == len(b):
        raise Exception('lengths of both inputs have to be the same')
    for i in range(len(a)):
        if not (a[i] == b[i]):
            print(i)
            print("python:", a[i])
            print("matlab:", b[i])
            x = False
    return x


def compare_arrays(a, b, round_val: int = -1, debug: bool = False) -> bool:
    """
    Compares an array of numbers ('a') against an array of arrays, where each sub-array contains only one number.
    For example, it compares [1, 2, 3] against [[1], [2], [3]].
    Args:
        a: Python array
        b: Matlab array
        round_val (int, optional): Determines whether the inputs should be rounded.
                                   Use -1 for no rounding (default), otherwise specify the number of decimals.
        debug (bool, optional): Enable debug mode. False by default.

    Returns: True if the input are equal, False otherwise.

    """
    if round_val > -1:
        a = np.round(a, round_val)
        for i in range(len(a)):
            b[i][0] = round(b[i][0], round_val)
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    x = True
    try:
        if a == b:
            return True
    except:
        logging.debug('Iterating through the entire array')
    if not len(a) == len(b):
        raise Exception('lengths of both inputs have to be the same')
    for i in range(len(a)):
        if not (a[i] == b[i][0]):
            print(i)
            print("python:", a[i])
            print("matlab:", b[i][0])
            x = False
    return x


if __name__ == "__main__":
    print("""You can't run this. This is a helper file""")
