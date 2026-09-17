# OSET ECG Feature Extraction Toolset

An open-source ECG feature extraction toolbox that incorporates an extensive list of clinically relevant ECG characteristics. The feature set can be applied to datasets with different ECG lengths (preferably 10 seconds or longer for valid statistical analysis), arbitrary lead configurations, and various sampling frequencies. This toolbox extracts a wide range of interpretable features from ECG signals for machine learning and statistical analysis. The extracted features are provided as a **.csv file**, along with detailed feature descriptions and corresponding units.

## Feature Categories:

The following categories include the ECG features provided by the toolbox. For detailed information about the feature definitions and extraction methods, please refer to our published papers listed below:

1. [**Heart Rate Variability and Heart Rate Metrics**](./ecg_hrv_features.m)
2. [**Beat signal-to-noise ratios**](./ecg_snr_features.m)
3. [**Singular Value Decomposition Metrics**](./ecg_svd_features.m)
4. [**Amplitude-to-Timing Ratios**](./ecg_amp_to_int_ratio_features.m)
5. [**Amplitude and Morphological Area Metrics**](./ecg_area_amp_features.m)  
6. [**Time Interval Measurements**](./ecg_time_intervals_features.m) 
7. [**Hjorth Descriptors**](./ecg_hjorth_features.m)
8. [**Average Beat Morphology**](./ecg_mean_phase_beat.m)
9. [**Cardiac Axis**](./ecg_angle_features.m)

## Requirements

- Input ECG signals should be in **WFDB** format and expressed in **millivolts (mv)** or **microvolts (μV)**.
- The code requires the **OSET package** for preprocessing, R peak detection, and fiducial point detection.

## Running the Code

The [**root folder**](./) folder contains all ECG feature-extraction functions you can use in your projects.
The [**test_script.m**](./run_tools/test_script_extract_features.m) provides an example of how to run all feature extraction functions. Example output files:

- [**lead-wise feature example**](../../../../datasets/sample-data/HR00001/leadwise_features/HR00001_features_leadwise.csv)
- [**multi-lead feature example**](../../../../datasets/sample-data/HR00001/multileads_features/HR00001_features_multileads.csv)

## Authors:
- [**Seyedeh Somayyeh Mousavi**](https://scholar.google.com/citations?user=gk99WMsAAAAJ&hl=en)
- [**Sajjad Karimi**](https://scholar.google.com/citations?user=nIUVskwAAAAJ&hl=en)
- [**Reza Sameni**](https://scholar.google.com/citations?user=MkoXtWwAAAAJ&hl=en)

## Citation
If you find this toolbox useful in your research or projects, please consider citing the publications associated with its development and application. Your citation helps acknowledge the work behind the toolbox and supports its continued development.

1. ***Estimating Blood Pressure from the Electrocardiogram: Findings of a Large-Scale Negative Results Study***,  
   Mousavi SS, Karimi S, Hassannia MS, Koscova Z, Rad AB, Albert D, Clifford GD, Sameni R. *Physiological Measurement*, Nov 2025, DOI: [10.1088/1361-6579/ae1926](https://doi.org/10.1088/1361-6579/ae1926).
2. ***Electrocardiogram-Based Machine Learning Model for Predicting Adverse Cardiac Outcomes in Adult Congenital Heart Disease***,
   Mousavi SS, Raskind-Hood CL, Haffner A, Robichaux C, Ivey LC, M Book W, Sameni R., Biomedical Signal Processing and Control, May 2026 (Under review).
3. ***ECG-Age for Cardiovascular Outcome Prediction: Evaluation in an Adult Congenital Heart Disease Cohort***,
   Mousavi SS, Raskind-Hood CL, Haffner A, M Book W, Sameni R., Computing in Cardiology (CinC) Conference, Sep 2026, [Link](https://cinc.org/2026/Program/accepted/287_Preprint.pdf).
4. ***The open-source electrophysiological toolbox (OSET)***, Sameni R. Version 4.0. 2006–2026, (https://github.com/alphanumericslab/OSET).
   

