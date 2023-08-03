# Python function that returns the beats representation of an ECG
# Inputs: ecg_bwr - the input ECG
#         mids - the middle points between R-peaks, which represents the boundaries between different beats
# Output: x_hhs - ECG beat representation


def ECGHeartBeatsList(ecg_bwr, mids):
    x_hbs = []
    for i in range(len(mids) - 1):
        x_hb = ecg_bwr[mids[i] : mids[i + 1]]
        x_hbs.append(x_hb)
    return x_hbs
