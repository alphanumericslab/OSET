# ECG QT interval estimator

Estimates QT interval and related parameters from single or multichannel ECG data.

## Introduction

The function [`qt_estimator_sum_of_gaussians.m`](./qt_estimator_sum_of_gaussians.m) estimates the QT-interval and related parameters from multichannel ECG data. It uses a sum-of-gaussian functions to model the ECG segments within a Maximum Likelihood (ML) and Bayesian framework. See reference below for the theoretical background and methods.

## Usage
```matlab
   fiducials = qt_estimator_sum_of_gaussians(ecg, fs, params);
```

## Parameters

The function accepts an optional parameter structure `params` with the following options:

- `w1`: Window length of the first median filter used in baseline removal. Default: 0.75
- `w2`: Window length of the second moving average filter used in baseline removal. Default: 0.9
- `f_hr`: Default heart rate factor for peak detection in Hz ( = BPM/60). Default: 1.2
- `amp_th`: Amplitude threshold for peak detection. Default: 0.6
- `polarity_detection_mode`: Polarity detection mode for peak detection. Default: 2
- `initial_q`: Initial values for Q wave segments of interest. Default: [-0.080; -0.020]
- `initial_t`: Initial values for T wave segments of interest. Default: [0.1; 0.5]
- `BayesianEstimation`: Perform Bayesian estimation or not (true/false). Saves processing load if only ML estimation is required.
- `snr`: Signal-to-noise ratio for Bayesian approach. Default: 100 (dB)
- `beta`: Constant coefficient for Q-T wave width. Default: 3. Can be a construct in this form:
   - `beta.q`: coefficient multiplied by the Q-wave gaussian STD to obtain the Q-wAve offset/onset (to be optimized per application)
   - `beta.t`: coefficient multiplied by the T-wave gaussian STD to obtain the T-wAve offset/onset (to be optimized per application)
   - `widest_qt`: widest QT-interval. QT estimates beyond this threshold are set to nan. Default: .6 seconds
   - `earliest_q_onset`: earliest Q-wave onset relative to the R-peak. Q-wave onsets and offsets are both set to nan, if the Q onset or offsets are estimated to be earlier. Default: 0.15 seconds
   - `latest_t_offset`: latest T-wave offset relative to the R-peak. T-wave onsets and offsets are both set to nan, if the T onset or offsets are estimated to be later. Default: 0.4 seconds
- `plot_results`: Whether to plot results. Default: false

## Outputs

The output structure `fiducials` contains various fiducial points and intervals estimated using both ML and Bayesian approaches. The `fiducials` structure includes the following fields, each of which is a vector of the same length as the number of detected R-peaks:

- `R_peaks`: Detected R peak positions for each channel.
- `QT_ML`: QT intervals estimated using the Maximum Likelihood approach for each channel.
- `Q_onsets_ML`: Q wave onset positions estimated using the Maximum Likelihood approach for each channel.
- `Q_offsets_ML`: Q wave offset positions estimated using the Maximum Likelihood approach for each channel.
- `T_onsets_ML`: T wave onset positions estimated using the Maximum Likelihood approach for each channel.
- `T_offsets_ML`: T wave offset positions estimated using the Maximum Likelihood approach for each channel.
- `QT_BYS`: QT intervals estimated using the Bayesian approach for each channel.
- `Q_onsets_BYS`: Q wave onset positions estimated using the Bayesian approach for each channel.
- `Q_offsets_BYS`: Q wave offset positions estimated using the Bayesian approach for each channel.
- `T_onsets_BYS`: T wave onset positions estimated using the Bayesian approach for each channel.
- `T_offsets_BYS`: T wave offset positions estimated using the Bayesian approach for each channel.


## Example:
```matlab
% Sample ECG data (replace with your actual ECG data)
load SampleECG1kHz2 data
data = data(:, 2:end)'; % remove the time column and make the data matrix in channels times samples format
fs = 1000;  % Sampling frequency

% Parameters (optional, modify as needed)
params.w1 = 0.75;
params.w2 = 0.9;
params.f_hr = 1.2;
params.amp_th = 0.6;
params.polarity_detection_mode = 2;
% default Q-wave/T-wave onsets and offsets with respect to R-peak position (only a rough estimate; will be optimized by the function):
params.initial_q = [-0.080; -0.020];  
params.initial_t = [0.1; 0.5];
params.BayesianEstimation = true; % run Bayesian estimator (in addition to ML)
params.snr = 10;  % Signal-to-noise ratio in dB (if known; used in the Bayesian estimator). Set to a high value such as 100 if unknown
params.beta.q = 3.1; % coefficient multiplied by the Q-wave gaussian STD to obtain the Q-wAve offset/onset (to be optimized per application)
params.beta.t = 2.9; % coefficient multiplied by the T-wave gaussian STD to obtain the T-wAve offset/onset (to be optimized per application)
params.widest_qt = 0.6; % widest QT
parsms.earliest_q_onset = 0.1; % earliest Q-wave onset
parsms.latest_t_offset = 0.35; % latest T-wave offset
params.plot_results = true; % plot the results

% Call qt_estimator_sum_of_gaussians to get the fiducial points
fiducials = qt_estimator_sum_of_gaussians(ecg_data, fs, params);

% Display the results
disp('QT Intervals (ML):');
disp(fiducials.QT_ML);

disp('QT Intervals (Bayesian):');
disp(fiducials.QT_BYS);
```

## Reference
Davood Fattahi and Reza Sameni. "Cramér-Rao Lower Bounds of Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions on Signal Processing 70 (2022): 3181-3192.

## Contributors
Davood Fattahi and Reza Sameni, 2023, [The Open-Source Electrophysiological Toolbox](https://github.com/alphanumericslab/OSET)