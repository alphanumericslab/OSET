# ECG Feature Extraction Pipeline

This project provides an ECG feature extraction pipeline that extracts features from ECG signals. The pipeline aims to report ECG signal features using a statistical approach (mean, median, and standard deviation for each feature). Additionally, there are other features to consider, including Singular Value Decomposition (SVD)-based features, that should be added to the 84 extracted features per lead.

## Overview

- The ECG signals should be in the **WFDB** format.
- The ECG signals should be in **millivolts (mv)** or **microvolts (μV)**.
- All signals should have the same number of leads and include noise filtering (notch filter) to remove power line noise.
- The output will be a **CSV file** that contains the extracted features, with units reported for each feature.

## Requirements

- **ECG signals in WFDB format**:  
  Ensure that the input ECG signals are provided in the WFDB format (.mat and .hea).

- **OSET Package**:  
  This project utilizes the OSET package for preprocessing, R peak detection, and fiducial point detection.  
  To use the package, follow these steps:  
  1. Clone the OSET repository
  2. Add the OSET package to your MATLAB path 

## Running the Code

To extract the features, use the [`extract_ecg_features_path_wfdb`](extract_ecg_features_path_wfdb.m) function. Here is a [`test script`](test_script.m) to run the function and this is an example of [`.csv output file`](ECG_features.csv). Below are the parameters you need to set before running the function:

### Parameters:
- **Path to the folder containing input ECG signals**  
  Define the path to the directory containing your ECG signals.

- **Path to save the output CSV file**  
  Specify the directory where the output CSV file will be saved.

- **Name of the output CSV file**  
  Define the name for the output CSV file that will contain the extracted features.

- **Lead List**  
  Define the list of ECG leads (e.g., `["I", "II", "III"]`).

- **Number of SVD Features to Compute per Lead**  
  Specify the number of Singular Value Decomposition (SVD) features to compute per lead.

- **Frequency for Notch Filtering (double)**  
  Define the notch filter frequency (in Hz), typically used to remove power line interference.

### Example Usage:

```matlab
% Set the parameters
input_path = 'path_to_input_folder';   % Path to folder containing ECG signals
output_path = 'path_to_output_folder'; % Path to save the output CSV file
output_filename = 'extracted_features.csv'; % Output CSV file name
lead_list = {'I', 'II', 'III'};       % List of ECG leads
n_svd_features = 10;                   % Number of SVD features to compute per lead
f_notch = 50;                          % Notch filter frequency (Hz)

% Call the feature extraction function
extract_ecg_features_one_record(input_path, output_path, output_filename, lead_list, n_svd_features, f_notch);
```

## ECG Features Extracted

1. [**SNR Features**](ecg_snr_features.m) (2):  
 The function calculates and reports the mean and median signal-to-noise ratio (SNR) based on noise levels for beats of lead.

2. [**SVD Features**](ecg_svd_features.m):  
  Singular Value Decomposition (SVD) features, with the number of features defined by the user.

3. [**HRV Features**](ecg_hrv_features.m) (7):  
  This function extracts Heart Rate Variability (HRV) features, including:
   - **SDNN** (Standard Deviation of NN intervals)  
   - **RMSSD** (Root Mean Square of Successive RR Interval Differences)  
   Additionally, the function calculates:
   - Median heart rate
   - Mean heart rate
   - Upper and lower 5% heart rate
   - Number of beats

4. [**Angles Features**](ecg_angles_features.m) (12):  
   The angles between R peaks and the P, S, Q, and T peaks are calculated, and the mean, median, and standard deviation of these angles are reported.

5. [**Amplitude and Area under the Curve Features**](ecg_area_amp_features.m) (32):  
   Features related to the ECG signal's morphology and time intervals:  
   - QRS complex, P-wave, and T-wave amplitudes  
   - Area under the curve for these components  
   - ST segment analysis,
   - Ratio of the amplitude of R to T and P waves.
   - Corrected QT interval, using the Fridericia and Bazett formulas.

6. [**Time Interval Features**](ecg_time_intervals_features.m) (29):  
   Features related to time intervals, reporting the mean, median, and standard deviation for all waves, segments, and intervals.

7. [**ECG Complexity Analysis**](ecg_complexity_features.m) (2):  
   Mobility and complexity calculations.


