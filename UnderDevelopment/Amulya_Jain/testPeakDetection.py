import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import PeakDetection

# Load the MATLAB .mat file
mat = scipy.io.loadmat('SampleECG1.mat')['data'][0]
fs = 1000
t = np.arange(len(mat)) / fs
f = 1

# Assuming PeakDetection is a function that returns peaks
# You need to implement PeakDetection separately or use an existing library for peak detection
peaks1 = PeakDetection.peak_detection(mat, f / fs)[0]
print(np.where(peaks1))
# Plotting
plt.figure()
plt.plot(t, mat, 'b', label='ECG')
plt.plot(t, peaks1 * mat, 'ro', label='ECG Peaks (max detection)')
plt.xlabel('time (sec.)')
plt.legend()
plt.grid()
plt.show()
