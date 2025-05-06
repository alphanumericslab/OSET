# ECG Feature Extraction Pipeline

This project is part of the OSET (Open-Source Electrophysiological Toolbox) and provides a comprehensive ECG feature extraction pipeline. The toolbox extracts a wide range of features from ECG signals that can be used for machine learning and statistical analysis applications.

## Overview

The pipeline processes ECG signals to extract features using a statistical approach (mean, median, and standard deviation for each feature). The extracted features can be used for:
- Machine learning model development
- Statistical analysis of ECG data
- Clinical research and diagnosis support

## Key Features
- Processes ECG signals in **WFDB** format (.mat or .dat files with accompanying .hea files)
- Handles signals in **millivolts (mV)** or **microvolts (μV)**
- Supports multi-lead ECG recordings
- Performs automatic noise filtering and preprocessing
- Exports features to **CSV files** with units for easy integration with ML frameworks
- Provides detailed feature descriptions

## Requirements

- **MATLAB** (R2019b or later recommended)
- **WFDB format ECG signals**:  
  Input ECG signals should be in WFDB format (.mat/.dat with .hea files)
- **OSET Package**:  
  The code requires the OSET package for preprocessing, R peak detection, and fiducial point detection

## Getting Started

### Setup
1. Clone the OSET repository or download this ECG feature extraction module
2. Add the OSET package to your MATLAB path:
   ```matlab
   addpath(genpath('/path/to/OSET/matlab/tools/ecg'));
   ```

### Running the Code

There are two main ways to use this feature extraction pipeline:

#### Option 1: Using `process_ecg_wfdb` (Recommended)

This function processes WFDB format ECG files and saves extracted features to CSV files.

```matlab
% Example usage:
input_wfdb_address = './sample-ecg/HR00001.mat';  % Path to ECG file
ecg_csv_file_name = './output/HR00001_features.csv';  % Output file
lead_names_target = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
window_dur = 10;    % Process in 10-second windows
start_time = 0;     % Start from beginning
stop_time = 0;      % Process until the end
flatten_flag = 0;   % Stack features by channel (0) or flatten (1)

process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time, flatten_flag);
```

For a complete demonstration with various ECG file formats, refer to [`demo_wfdb_process_ecg.m`](demo_wfdb_process_ecg.m).

#### Option 2: Using `ecg_feature_extraction` (For advanced users)

This approach gives you more control over the preprocessing and feature extraction pipeline.

```matlab
% Load and preprocess your data
[ecg_data, fs, lead_names] = load_your_ecg_data();  % Replace with your loading function

% Extract features
flatten_flag = false;  % Features stacked by channel
[ecg_features_vector, ecg_feature_info, ecg_fiducial_position] = ecg_feature_extraction(ecg_data, fs, lead_names, [], [], [], [], flatten_flag);

% Save features to CSV
feature_handle_csv('output_features.csv', ecg_features_vector, ecg_feature_info.names, ecg_feature_info.units, ecg_feature_info.description, fs, window_time_info, lead_names);
```

For a complete example of preprocessing and feature extraction, refer to [`sample_run_ecg_feature_extraction.m`](sample_run_ecg_feature_extraction.m).

### Important Parameters

- **input_wfdb_address**: Path to the WFDB format ECG file
- **ecg_csv_file_name**: Path to save the output CSV file
- **lead_names_target**: List of ECG leads to process
- **window_dur**: Duration of windows for processing (in seconds)
- **flatten_flag**: 
  - `0`: Features are organized by channel [Channels × Features]
  - `1`: Features are flattened across channels [1 × (Channels × Features)]

## Feature Categories

The pipeline extracts the following categories of features:

1. **SNR Features** (2): Signal-to-noise ratio measurements
2. **SVD Features**: Singular Value Decomposition features (user-defined count)
3. **HRV Features** (7): Heart Rate Variability metrics (SDNN, RMSSD, etc.)
4. **Angles Features** (12): Angular relationships between different ECG waves
5. **Amplitude and Area Features** (32): Wave amplitudes and areas under the curve
6. **Time Interval Features** (29): Durations of waves, segments, and intervals
7. **Complexity Features** (2): Mobility and complexity measurements
8. **Morphology Features**: Sampled values from the average beat (user-defined count)

## Preparing Features for Machine Learning

The features exported to CSV can be directly used for machine learning applications:

1. **Feature Selection**: The CSV includes feature names and descriptions to help you select relevant features
2. **Feature Engineering**: You can combine or transform features as needed
3. **Model Training**: Import the CSV into your preferred ML framework

```python
# Example of using the features in Python (after extraction)
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

# Load features
features = pd.read_csv('output_features.csv')
X = features.iloc[:, 2:]  # Assuming first two columns are metadata
y = your_labels  # Your classification labels

# Split and train
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
model = RandomForestClassifier()
model.fit(X_train, y_train)
```

## Advanced Usage

### Processing Large Datasets

For large datasets, we recommend:
1. Processing files in batches
2. Using the windowing approach with appropriate `window_dur`
3. Parallel processing multiple files when possible

### Custom Feature Subsets

You can extract specific feature subsets by modifying the feature extraction function calls or by filtering the output CSV.

## Troubleshooting

- **Missing WFDB Header**: Ensure both .dat/.mat and .hea files are present
- **Lead Naming**: Lead names are case-sensitive and must match those in the header
- **Memory Issues**: Reduce window size for large recordings
- **NaN Values**: Check signal quality or try different preprocessing parameters

## Further Information

For more details on specific feature extraction algorithms, refer to the documentation in the individual feature extraction scripts.

## Authors & Contact

### Contributors
- **Sajjad Karimi**
  - Location: Emory University, Georgia, USA
  - Email: sajjadkarimi91@gmail.com
  
- **Seyedeh Somayyeh Mousavi**
  - Location: Emory University, Georgia, USA
  - Email: bmemousavi@gmail.com

### Project Lead
- **Reza Sameni**
  - Email: reza.sameni@gmail.com

For questions, issues, or collaborations related to the OSET toolbox, please contact the authors.


