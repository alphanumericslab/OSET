# Python function that is reording the R-peaks detected on the signal in order to assure "full" beats.
# Input: peaks - the R-peaks positions
#     : mids - the mid points positions
# Output: pks, the reordered list of peaks.


def RPeakReordering(peaks, mids):
    if peaks[0] < mids[0] and peaks[-1] > mids[-1]:
        pks = peaks[1:-1]
        # print("case 1")
    elif peaks[0] < mids[0]:
        pks = peaks[1:]
        # print("case 2")
    elif peaks[-1] > mids[-1]:
        pks = peaks[:-1]
        # print("case 3")
    else:
        pks = peaks
    return pks
