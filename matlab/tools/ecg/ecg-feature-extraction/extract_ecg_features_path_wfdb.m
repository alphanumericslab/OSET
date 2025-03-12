function extract_ecg_features_path_wfdb(input_path, output_path, file_name, lead_list, n_eigenvalues, f_notch)
% extract_ecg_features_path_wfdb(input_path, output_path, file_name, lead_list, n_eigenvalues, f_notch)
% Description: Extracts features from ECG signals in WFDB format and saves the results to a .csv file.
% 
% Inputs:
%   input_path     - Path to the folder containing input (ECG signals in WFDB format) (string)
%   output_path   - Path to save the output .csv file (string)
%   file_name      - Name of the output .csv file (string)
%   lead_list        - Cell array specifying the ECG leads to analyze (e.g., {'I', 'II', 'V1', ...})
%   n_eigenvalues  - Number of eigenvalues (SVD features) to compute per lead (integer)
%   f_notch        - Frequency for notch filtering (double)
%
% Outputs:
%   None. The extracted features are saved directly to the specified .csv file.
%
% Notes:
%   - Requires WFDB files in the specified input directory.
%   - The output file includes a structured header with units and features.
%
% Date: Jan 10, 2025
% Author: Seyedeh Somayyeh Mousavi
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com

%% constant  value
convert_mv_um = 1000;

%% Ensure the output directory exists
if ~isfolder(output_path)
    error('Output path does not exist: %s', output_path);
end

% Open the CSV file for writing
csv_file = fopen(fullfile(output_path, file_name), 'w');
if csv_file == -1
    error('Failed to open file: %s', fullfile(output_path, file_name));
end

%% Define the static header
fprintf(csv_file, 'filename, ECG_length, N_lead, Unit, Gain');

% Define ECG-related features
ecg_features = {
    'SNR_median', 'SNR_mean', 'N_beats', 'RMSSD', 'SDNN', 'HR_median', 'HR_mean', 'HR_upper_5', 'HR_lower_5', ...
    'mean_pr_interval_TI', 'std_pr_interval_TI', 'median_pr_interval_TI', 'mean_qt_interval_TI', 'std_qt_interval_TI', 'median_qt_interval_TI', 'QTc_b','QTc_f' ...
    'mean_st_Se_TI', 'std_st_Se_TI', 'median_st_Se_TI', 'mean_pr_Se_TI', 'std_pr_Se_TI', 'median_pr_Se_TI', ...
    'mean_pr_peaks_TI', 'std_pr_peaks_TI', 'median_pr_peaks_TI', 'mean_qr_peaks_TI', 'std_qr_peaks_TI', 'median_qr_peaks_TI', ...   
    'mean_qrs_complex_TI', 'std_qrs_complex_TI', 'median_qrs_complex_TI', ...
    'mean_rs_peaks_TI', 'std_rs_peaks_TI', 'median_rs_peaks_TI', 'mean_rt_peaks_TI', 'std_rt_peaks_TI', 'median_rt_peaks_TI', ...       
    'mean_qrs_A', 'std_qrs_A', 'median_qrs_A', 'mean_qrs_area', 'std_qrs_area', 'median_qrs_area', 'mean_qrs_abs_area', 'std_qrs_abs_area', 'median_qrs_abs_area', ...
    'mean_t_A', 'std_t_A', 'median_t_A', 'ratio_qrs_t_A', 'mean_t_area', 'std_t_area', 'median_t_area',  'mean_t_abs_area', 'std_t_abs_area', 'median_t_abs_area', ...
    'mean_p_A', 'std_p_A', 'median_p_A', 'ratio_qrs_p_A', 'mean_p_area', 'std_p_area', 'median_p_area', 'mean_p_abs_area', 'std_p_abs_area', 'median_p_abs_area', ...
    'mean_st_Se_A', 'std_st_Se_A', 'median_st_Se_A', ...
    'mean_pr_angles', 'std_pr_angles', 'median_pr_angles', ...
    'mean_qr_angles', 'std_qr_angles', 'median_qr_angles', ...
    'mean_rs_angles', 'std_rs_angles', 'median_rs_angles', ...
    'mean_rt_angles', 'std_rt_angles', 'median_rt_angles' , ...
    'mobility', 'complexity'};

% Define corresponding units for ECG features
ecg_units = {
    'db', 'db', 'scaler', 'ms', 'ms', 'bpm', 'bpm', 'bpm', 'bpm',...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms','ms', 'ms', ...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms', ...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms', ...
    'ms', 'ms', 'ms', ...
    'ms', 'ms', 'ms', 'ms', 'ms', 'ms', ...
    'mv', 'mv', 'mv', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', ...
    'mv', 'mv', 'mv', 'scaler', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms',...
    'mv', 'mv', 'mv', 'scaler', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms', 'mv*ms',...
    'mv', 'mv', 'mv',...
    'deg', 'deg', 'deg', 'deg', 'deg', 'deg', ...
    'deg', 'deg', 'deg', 'deg', 'deg', 'deg', ...
    'scaler', 'scaler'};

% Define formatting rules based on the units
ecg_strings = cell(size(ecg_units));
for i = 1:length(ecg_units)
    switch ecg_units{i}
        case 'db'
            ecg_strings{i} = '%.3f';
        case 'bpm'
            ecg_strings{i} = '%.1f';
        case 'ms'
            ecg_strings{i} = '%.1f';
        case {'mv', 'scaler'}
            ecg_strings{i} = '%.3f';
        case 'mv*ms'
            ecg_strings{i} = '%.4f';
        case 'deg'
            ecg_strings{i} = '%.2f';
    end
end

% Update formatting for N_beats
index = find(strcmp(ecg_features, 'N_beats'));
if ~isempty(index)
    ecg_strings{index} = '%d'; % Use integer format for N_beats
else
    error('Feature "N_beats" not found in ecg_features');
end

% Append formatting for eigenvalues
ecg_strings = [ecg_strings, repmat({'%.2f'}, 1, n_eigenvalues)];

% Write feature names for each lead
for lead = 1:length(lead_list)
    for i = 1:length(ecg_features)
        fprintf(csv_file, ',%s_lead_%s', ecg_features{i}, lead_list{lead});
    end
    for i = 1:n_eigenvalues
        fprintf(csv_file, ',S_%d_lead_%s', i, lead_list{lead});
    end
end

% Add a newline to separate the header from the data
fprintf(csv_file, '\n');

% Write unit headers
fprintf(csv_file, '-, s, scaler, -, scaler');
for lead = 1:length(lead_list)
    for i = 1:length(ecg_units)
        fprintf(csv_file, ',%s', ecg_units{i});
    end
    for i = 1:n_eigenvalues
        fprintf(csv_file, ',scaler');
    end
end
fprintf(csv_file, '\n');

% Repeat ecg_strings for all leads
ecg_strings = repmat(ecg_strings, 1, length(lead_list));
%% ==================================================================================================
% Get a list of all .mat files in the folder and subfolders (ECG signals)
file_list = dir(fullfile(input_path, '**', '*.mat'));
%% ==================================================================================================
fprintf('Loop through each file, process the data, and extract features.\n')
fprintf('Number of ECG files in the folder: %d\n', length(file_list))

n_features_ecg = size(ecg_features, 2);
n_features = (n_features_ecg + n_eigenvalues )* length(lead_list);

% % Loop through each file, process the data and extract features
for i = 1:length(file_list)

      % Display only when 'i' is a multiple of 200
       if mod(i, 200) == 0
          disp(i)  
       end

       % Use fileparts to get the base name and path
       [~, base_name, ~] = fileparts(file_list(i).name);
       data_path = file_list(i).folder;
       full_file_path = fullfile(data_path, base_name);

       % load the signal
       [tm,data,fs ,siginfo]= rdmat(full_file_path);

       sig_len = length(tm)/fs;
       unit = {siginfo.Units};
       gain = {siginfo.Gain};
       challenge_str = {siginfo.Description};
 
       % Check if data has the right dimensions
       if size(data, 1) > size(data, 2)             % Check if rows > columns
              data = data';                                  % Transpose the data
       end

       n_channels = size(data, 1);

       for h = 1:n_channels

                        if ismember(lower(unit{h}), {'mv'})                 % Check if unit_PART is 'mV', 'mv', 'MV', or 'Mv'
                                convert_ratio = convert_mv_um;    
                        elseif ismember(lower(unit{h}), {'uv'})            % Check if unit_PART is 'uV', 'uv', 'UV', or 'Uv'            
                                convert_ratio = 1;
                        else           
                        % Handle other cases if needed
                        disp('Unit is not: mv and uv');
                        fprintf("It seems there is an issue with the gain part of the signal: %s, channel: %d \n", base_name, n_channels );
                        continue;
                        end
                        data(h, :) = convert_ratio *  data(h, :);
        end  
        % Check if all elements of unit are the same
        if all(strcmp(unit, unit{1}))
               commonUnit = unit{1};
        else
               commonUnit = NaN; 
        end

        % Check if all elements of gain are the same
        if all(cellfun(@(x) isequal(x, gain{1}), gain))
               commonGain = gain{1};
        else
               commonGain = NaN; 
        end
        
        % Reorder data if the lead order does not match
        if ~isequal(lower(challenge_str), lower(lead_list))
            
            % Find the indices in challenge_str that match n_lead_list
            [~, new_order] = ismember(lead_list, challenge_str);
    
            % Check for invalid indices (leads not found)
            if any(new_order == 0)
                fprintf("It seems there is an issue with the channel order of the signal: %s \n", base_name2);
                continue;  % Skip the current iteration and move to the next file
            end
    
            % Reorder the data according to the new order
            data = data(new_order, :);
        end

        all_features = extract_ecg_features_one_record(base_name, data, fs, n_eigenvalues, n_features, f_notch);
        
        % Check if all values in all_features are NaN
        if all(isnan(all_features))
                fprintf('All features are NaN. Skipping save operation.\n');
        else
                % Write the filename and data length to the CSV file
                fprintf(csv_file, '%s, %d, %d, %s, %d,', base_name, sig_len, n_channels, commonUnit, commonGain);

                % Write all concatenated features to the CSV file
                for k = 1:length(all_features)
                        fprintf(csv_file, ecg_strings{k}, all_features(k)); 
                         if k < length(all_features)
                                fprintf(csv_file, ',');  % Add a comma between values
                        end

                end
                % Move to the next line after writing the row
                fprintf(csv_file, '\n');
         end
end
fprintf('Extracting feature is finished.\n')
        
% Close the CSV file
fclose(csv_file);
fprintf('%s file is saved.\n', file_name);

