function extract_ecg_features( input_path, output_path, opts)

% function extract_ecg_features(input_path, output_path, opts)
% 
% Inputs:
%   input_path (string or char): Path to the folder containing input ECG signals (WFDB format)
%   output_path (string or char): Path to save the extracted feature output (.csv file)
%
% Outputs:
%   None. The extracted features are saved directly as a .csv file in the output_path.
%
% Notes:
%   Requires WFDB files in the input directory.
%
% Authors: 
%   Seyedeh Somayyeh Mousavi
%   Reza Sameni
%
% Emory University, Georgia, USA
% bmemousavi@gmail.com
% AUG, 2026

% Features
%   1. HRV features (12) (with N_beats and ECG_length)
%   2. SNR features (2)
%   3. Mean ECG beat (n_mean_ECG_beats)
%   4. Features related to time interval (29)
%   5. Features related to amplitude and area under curve (32)
%   6. Features related to amplitude to time interval ratios (12)
%   7. Complexity analysis features (3)
%   8. SVD features (n_svd)
%   9. ECG angle features (18)

arguments

    input_path  (1,:) char {mustBeFolder}
    output_path (1,:) char
    opts.lead_list (1,:) cell = ...
                        {"I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"};
    opts.N_svd (1,1) double {mustBeInteger, mustBePositive} = 5
    opts.N_mean_ECG_beats (1,1) double {mustBeInteger, mustBePositive} = 100
    opts.flag_baseline_alignment (1,1) double {mustBeMember(opts.flag_baseline_alignment,[0,1])} = 0
    opts.f_notch (1,1) double {mustBePositive} = 60
    opts.low_pass_fre_leadwise (1,1) double {mustBePositive} = 0.1
    opts.high_pass_fre_leadwise (1,1) double {mustBePositive} = 100
    opts.low_pass_fre_multileads (1,1) double {mustBePositive} = 0.5
    opts.high_pass_fre_multileads (1,1) double {mustBePositive} = 40
    opts.median_window_coff (1,1) double {mustBePositive} = 0.3
    opts.mean_window_coff (1,1) double {mustBePositive} = 0.15
    opts.Q_factor (1,1) double {mustBePositive} = 40
    opts.Pre_R_peak_segment (1,1) double {mustBeInRange(opts.Pre_R_peak_segment,0,1,"exclusive")} = 0.3
    opts.Post_R_peak_segment (1,1) double {mustBeInRange(opts.Post_R_peak_segment,0,1,"exclusive")} = 0.7
    opts.increasing_ECG_beat_lenght (1,1) double {mustBePositive} = 1.2
    opts.flag_post_processing (1,1) double {mustBeMember(opts.flag_post_processing,[0,1])} = 1

end
rng(42)

%%  Check the values
if ~exist(output_path, 'dir')
     fprintf('Output path "%s" does not exist. Creating...\n', output_path);
     mkdir(output_path);
end

if (opts.Pre_R_peak_segment + opts.Post_R_peak_segment) ~= 1
    error('Pre_R_peak_segment and Post_R_peak_segment must sum to 1.');
end

%% Create output paths
leadwise_output_path    = fullfile(output_path, 'leadwise_features');
multileads_output_path  = fullfile(output_path, 'multileads_features');

if ~exist(leadwise_output_path, 'dir')
     mkdir(leadwise_output_path);
end

if ~exist(multileads_output_path, 'dir')
     mkdir(multileads_output_path);
end
%% ====================================================================
% Get a list of all .mat files in the folder and subfolders (ECG signals)
file_list = dir(fullfile(input_path, '**', '*.mat'));
%% ====================================================================
fprintf('Loop through each file, process the data, and extract features.\n')
fprintf('Number of ECG files in the folder: %d\n', length(file_list));

% Loop through each file
parfor i = 1:length(file_list)

    % Display progress every 200 files
    if mod(i, 200) == 0
           fprintf('Processing file %d of %d\n', i, length(file_list));
    end

    try 
           % Load record 
           [data, fs, base_name] = load_and_preprocess_ecg(file_list(i), opts.lead_list);
           
           
           %  Preprocess data
           data_filtered_for_leadwise_feature = ...
           preprocess_ecg_leads_for_leadwise_features(data,fs, opts.f_notch, opts.Q_factor, ...
                                                      opts.low_pass_fre_leadwise, opts.high_pass_fre_leadwise, ...
                                                      opts.median_window_coff, opts.mean_window_coff );  
           
           % Extract leadwise features
           extract_ecg_features_leadwise( ...
                        base_name, data_filtered_for_leadwise_feature, fs, ...
                        opts.increasing_ECG_beat_lenght, ...
                        opts.flag_post_processing, ...
                        opts.Pre_R_peak_segment, ...
                        opts.Post_R_peak_segment,...
                        opts.N_mean_ECG_beats, ...
                        opts.flag_baseline_alignment,  ...
                        opts.N_svd, ...
                        opts.lead_list, ...
                        leadwise_output_path);

           % Extract multi lead features    
           extract_ecg_features_multileads( ...
                        base_name, data, fs, ...
                        opts.low_pass_fre_multileads, ...
                        opts.high_pass_fre_multileads, ...
                        opts.median_window_coff, ...
                        opts.mean_window_coff, ...
                        opts.lead_list, ...
                        opts.flag_post_processing, ...
                        multileads_output_path);
    catch 

           fprintf("Error processing file %s: %s\n", file_list(i).name);
        
    end

end
