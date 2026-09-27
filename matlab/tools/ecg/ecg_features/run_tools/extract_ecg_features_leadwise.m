function [Table_csv] = extract_ecg_features_leadwise(signal_name, data, fs , ...
                                                                        increasing_ECG_beat_length, ...
                                                                        flag_post_processing, ...
                                                                        pre_R_peak_segment, post_R_peak_segment, ...
                                                                        n_mean_ECG_beats, flag_baseline_alignment, ...
                                                                        n_svd,  lead_list, output_path)

    % function [Table_csv] = extract_ecg_features_leadwise(signal_name, data, fs , ...
    %                                                         increasing_ECG_beat_lenght, ...
    %                                                         flag_post_processing, ...
    %                                                         pre_R_peak_segment, post_R_peak_segment, ...
    %                                                         n_mean_ECG_beats, flag_baseline_alignment, ...
    %                                                         n_svd,  lead_list, output_path)                                                  
    % Extract of all features from one record of ECG signal
    %
    % Inputs:
    %   signal_name (string): Name of the ECG signal
    %   data (1D array, microvolts): Preprocessed ECG signal
    %   fs (scalar, Hz): Sampling frequency of the ECG signal
    %   increasing_ECG_beat_lenght - (float) (Required for ecg_snr_features function)
    %   flag_post_processing (binary)(Required for the  function)
    %   pre_R_peak_segment (float): Portion used before the R-peak
    %   post_R_peak_segment (float): Portion used after the R-peak
    %   n_mean_ECG_beats (scalar): Number of samples to consider for morphology-based features    
    %   flag_baseline_alignment (binary)(Required for ecg_mean_phase_beat function)
    %   n_svd (scalar): Number of eigenvalues to consider for SVD-based features
    %   lead_list (cell array of strings): Names and order of ECG leads
    %   output_path (string): Full path to save the extracted features as a `.csv` file   
    %
    % DEPENDENCIES:
    %  1. peak_det_likelihood_long_recs function
    %  2. ecg_delineate_lsim function
    %
    % OUTPUT:
    % Table_csv: A csv file that includes the following features:
    %   1. HRV features (12) (with N_beats and ECG_length)
    %   2. SNR features (2)
    %   3. Mean ECG beat (n_mean_ECG_beats)
    %   4. Features related to time interval (29)
    %   5. Features related to amplitude and area under curve (32)
    %   6. Features related to amplitude to time interval ratios (12)
    %   7. Complexity analysis features (3)
    %   8. SVD features (n_svd)
    %
    % Author: 
    %   Seyedeh Somayyeh Mousavi
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA
    
    %% Constant values
    n_features = 12 + 2 + 29 + 32 + 12 + 3;
    [n_channels, ECG_points] = size(data);
    n_features =  n_features + n_svd + n_mean_ECG_beats;
    all_features_vector = nan(n_channels, n_features);

    %% Feature Extraction
    try
         seg_len_time = min(10, ECG_points/fs);    
         for j = 1:n_channels

                         data_channel = data(j, :);                          
                         % Run R peak detection
                         if seg_len_time <= 10
                                    [R_locs_indexes, R_peaks_indexes] = peak_det_likelihood(data_channel, fs);
                         else  
                                    [R_locs_indexes, R_peaks_indexes] = peak_det_likelihood_long_recs(data_channel, fs , seg_len_time);
                         end

                         % Run ECG fiducial points detector    
                         position = ecg_delineate_lsim(data_channel, fs, R_peaks_indexes, flag_post_processing);
                         % Extract features
                         [hrv_features, hrv_feature_info]  = ecg_hrv_features(data_channel, R_peaks_indexes, fs);
                         [snr_features, snr_feature_info, ~]  = ecg_snr_features(data_channel, R_peaks_indexes, increasing_ECG_beat_length, pre_R_peak_segment, post_R_peak_segment);
                         [ti_features, ti_feature_info]  = ecg_time_intervals_features(position, fs);
                         [amps_areas_features, amps_areas_feature_info]  = ecg_area_amp_features(data_channel, position, fs);
                         [amp_to_int_ratio_features, amp_to_int_ratio_feature_info]  = ecg_amp_to_int_ratio_features(data_channel, position, fs);
                         [hjorth_features,  hjorth_feature_info]  = ecg_hjorth_features(data_channel, fs);
                         [mean_phase_beat_features, mean_beat_feature_info] = ecg_mean_phase_beat(data_channel, R_locs_indexes, n_mean_ECG_beats, flag_baseline_alignment);
                         [svd_features, svd_feature_info]  = ecg_svd_features(data_channel, R_peaks_indexes, n_svd);

                         % Concatenate extracted features into a single vector
                         all_features = [hrv_features, snr_features,  ti_features, amps_areas_features, ...
                                 amp_to_int_ratio_features, hjorth_features , mean_phase_beat_features, svd_features];

                         % Store in output matrix
                         all_features_vector(j, :) = all_features;
          end

          % Check if all values in all_features are NaN
         if all(isnan(all_features_vector))
             fprintf("All features are NaN. Skipping save operation: %s\n", signal_name);
         else

             % Define info structs
             feature_infos = { ...
                                        hrv_feature_info, ...
                                        snr_feature_info, ...
                                        ti_feature_info, ...
                                        amps_areas_feature_info, ...
                                        amp_to_int_ratio_feature_info, ...
                                        hjorth_feature_info, ...
                                        mean_beat_feature_info, ...
                                        svd_feature_info ...
                                        };
        
             % Initialize combined metadata
             all_feature_names = {};
             all_features_units = {};
             all_feature_description = {};
        
             % Concatenate metadata from each struct
             for i = 1:numel(feature_infos)
                     fi = feature_infos{i};
                     all_feature_names = [all_feature_names, fi.names];
                     all_features_units = [all_features_units, fi.units];
                     all_feature_description = [all_feature_description, fi.description];
             end
                        
             % rows
             feature_names = all_feature_names(:);  
             % converts string-in-cell -> char-in-cell
             feature_names = cellstr(feature_names(:)');   
             feature_units = all_features_units(:);           
             feature_description = all_feature_description(:);
             feature_values = num2cell(all_features_vector); 
        
             % columns 
             feature_data = [
                                feature_description(:)'; 
                                feature_units(:)'; 
                                feature_values
                                ];
        
             % Row names (vertical)
             row_names = [{'Descriptions'; 'Units'}; lead_list(:)];
        
             % Table
             Table_csv = cell2table(feature_data, 'VariableNames', cellstr(feature_names(:)'));
        
             % Construct the full file path using output_path and signal_name
             csv_filename = fullfile(output_path, [signal_name '_features_leadwise.csv']);
        
             % Convert the table to a format that can be saved as CSV
             Table_csv = [table(row_names, 'VariableNames', {'Features'}), Table_csv];
        
             % Write the table to CSV
             writetable(Table_csv, csv_filename);

         end 

     catch 

         % Error handling
         fprintf("\nERROR while processing\n");
         fprintf("Signal : %s\n", signal_name);
         fprintf("Channel: %d\n", j);

    end         
end

