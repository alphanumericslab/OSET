function [feature_vec, feature_info] = ecg_angle_features(data, fs, ...
                                                                   low_pass_fre, high_pass_fre, ...
                                                                   median_window_coff, mean_window_coff, lead_list, ...
                                                                   flag_post_processing)

    %  [feature_vec, feature_info] = ecg_angle_features(data, fs, 
    %                                                            low_pass_fre, high_pass_fre, 
    %                                                            median_window_coff, mean_window_coff, lead_list, 
    %                                                            flag_post_processing)  
    %
    %  Inputs:
    %   data (channels x samples double): raw ECG signal
    %   fs (scalar): Sampling frequency of the ECG signal 
    %   low_pass_fre (scalar): low-pass cut-off frequency (Hz)
    %   high_pass_fre (scalar): high-pass cut-off frequency (Hz)
    %   median_window_coff (scalar): coefficient for the moving median window
    %   mean_window_coff (scalar): coefficient for the moving mean window
    %   lead_list (cell array of strings): Names and order of ECG leads
    %   flag_post_processing (binary)(Required for the ecg_delineate_lsim function)
    %
    % Outputs:
    %   feature_vec: A vector contains features related to ECG angles
    %   feature_info: A structure contains feature descriptions, names and units.   
    % 
    % Ref: 
    %   1- Mueller-Leisse, Johanna, et al. "Determining the QRS axis: visual estimation is equal
    %   to calculation." Herzschrittmachertherapie+ Elektrophysiologie 36.1 (2025): 70-74.
    %   2- Novosel, Dragutin, Georg Noll, and Thomas F. Lüscher. 
    %   "Corrected formula for the calculation of the electrical heart axis." 
    %   Croatian medical journal 40 (1999): 77-79.
    %
    % Dependencies:
    %   1. `mean_median_ecg_amp` function
    %   2. `calculate_ecg_angle` function
    %   3. `preprocess_for_ecg_angle_features` function
    %   4. `peak_det_likelihood_long_recs` function
    %   5. `peak_det_likelihood` function
    %   6. `ecg_delineate_lsim` function
    %  
    % Authors:
    %   Seyedeh Somayyeh Mousavi 
    %   Reza Sameni
    %
    % bmemousavi@gmail.com
    % Aug 2026
    % Emory University, Georgia, USA
    %%
    desired_leads = {"I", "II", "aVF"};
    lead_list_flat = strtrim(string(lead_list));   
    try 

            [~, desired_list] = ismember(lower(string(desired_leads)), lower(string(lead_list)));
                
            %% Calculate the mean and median amplitude
            if all(desired_list ~= 0)

                    % Dictionaries for storing 
                    mean_amplitude = dictionary();
                    median_amplitude = dictionary();

                    %% Data preprocessing
                    lead_idx_in_original = zeros(1, numel(desired_leads));
                    for k = 1:numel(desired_leads)
                                idx_k = find(strcmpi(lead_list_flat, desired_leads{k}));
                                lead_idx_in_original(k) = idx_k;
                    end
                    data_subset = data(lead_idx_in_original, :);
                    data_filtered = preprocess_for_ecg_angle_features(data_subset, ...
                           fs, low_pass_fre, high_pass_fre, median_window_coff, mean_window_coff); 

                    % R-peaks of the ECG lead (I)
                    data_filtered_I = data_filtered(1, :);

                    % Run R peak detection
                    seg_len_time = length(data_filtered_I)/fs;
                    if seg_len_time <= 10
                            [~, R_peaks_indexes] = peak_det_likelihood(data_filtered_I, fs);
                    else  
                            [~, R_peaks_indexes] = peak_det_likelihood_long_recs(data_filtered_I, fs, seg_len_time);
                    end

                    for j = 1:length(desired_leads)
                                          
                              lead_name = desired_leads{j};
                              data_channel = data_filtered(j, :);

                              % Run ECG fiducial points detector    
                              position = ecg_delineate_lsim(data_channel, fs, R_peaks_indexes, flag_post_processing);
                              Length_data = length(data_channel);
                              [mean_amp, median_amp] = mean_median_ecg_amp(data_channel, position, Length_data);

                              % Store results using lead name
                              mean_amplitude(lead_name)  = mean_amp;
                              median_amplitude(lead_name) = median_amp;

                    end
                
                    %% Formula: 
                    % 1 - atan2(aVF,I) (degree)   
                    % 2 - atan2(2*aVF, (sqrt(3)*I)) (degree) 
                    % 3 - atan2(2*II-I, sqrt(3)*I)(degree) 
        
                    corrected_ratio_2 = 2/sqrt(3);
                    corrected_ratio_31 = 2;
                    corrected_ratio_32 = sqrt(3);
                    waves = ["P", "R", "T"];

                    % Mean Amplitude    
                    for w = 1:length(waves)
        
                            wave = waves(w);
                            % Form 1
                            axis_mean.(wave).Form_1 = calculate_ecg_angle(mean_amplitude("aVF").(wave), mean_amplitude("I").(wave));
                            % Form 2
                            axis_mean.(wave).Form_2 = calculate_ecg_angle(corrected_ratio_2 * mean_amplitude("aVF").(wave), mean_amplitude("I").(wave));
                            % Form 3
                            axis_mean.(wave).Form_3 = calculate_ecg_angle(...
                                (corrected_ratio_31 * mean_amplitude("II").(wave) - mean_amplitude("I").(wave)), ...
                                 corrected_ratio_32 * mean_amplitude("I").(wave));
                    end
        
                    % Median Amplitude 
                    for w = 1:length(waves)
                            
                            wave = waves(w);
                            % Form 1
                            axis_median.(wave).Form_1 = calculate_ecg_angle(median_amplitude("aVF").(wave), median_amplitude("I").(wave));
                            % Form 2
                            axis_median.(wave).Form_2 = calculate_ecg_angle(corrected_ratio_2 * median_amplitude("aVF").(wave), median_amplitude("I").(wave));
                            % Form 3
                            axis_median.(wave).Form_3 = calculate_ecg_angle(...
                                (corrected_ratio_31 * median_amplitude("II").(wave) - median_amplitude("I").(wave)), ...
                                 corrected_ratio_32 * median_amplitude("I").(wave));
                    end

                    %% feature_vec  
                    feature_vec = [ ...
                            axis_mean.T.Form_1, axis_mean.T.Form_2, axis_mean.T.Form_3, ...
                            axis_mean.P.Form_1, axis_mean.P.Form_2, axis_mean.P.Form_3, ...
                            axis_mean.R.Form_1, axis_mean.R.Form_2, axis_mean.R.Form_3, ...
                            axis_median.T.Form_1, axis_median.T.Form_2, axis_median.T.Form_3, ...
                            axis_median.P.Form_1, axis_median.P.Form_2, axis_median.P.Form_3, ...
                            axis_median.R.Form_1, axis_median.R.Form_2, axis_median.R.Form_3 ...
                            ];
                    % Convert Inf to NaN
                    feature_vec(~isfinite(feature_vec)) = NaN;
                    feature_vec = round(feature_vec); 
            
            end        
    catch 

          fprintf('Error in ecg_angle_features: %s\n', ME.message);
          fprintf('  Identifier: %s\n', ME.identifier);
          for k = 1:numel(ME.stack)
                     fprintf('  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
          end
          feature_vec = NaN(1, 18);    
    end 

    % Define feature info
    feature_info.names = { ...
            "Taxis_mean.Form_1",    "Taxis_mean.Form_2",    "Taxis_mean.Form_3", ...
            "Paxis_mean.Form_1",    "Paxis_mean.Form_2",   "Paxis_mean.Form_3", ...
            "Raxis_mean.Form_1",    "Raxis_mean.Form_2",   "Raxis_mean.Form_3", ...
            "Taxis_median.Form_1",  "Taxis_median.Form_2", "Taxis_median.Form_3", ...
            "Paxis_median.Form_1", "Paxis_median.Form_2", "Paxis_median.Form_3", ...
            "Raxis_median.Form_1", "Raxis_median.Form_2", "Raxis_median.Form_3" ...
            };
    feature_info.units = repmat({"degree"}, 1, length(feature_vec));

    feature_info.description = { ...
    "T-axis using mean amplitude: atan2(aVF, I)", ...
    "T-axis using mean amplitude: atan2(2*aVF, sqrt(3)*I)", ...
    "T-axis using mean amplitude: atan2(2*II-I, sqrt(3)*I)", ...
    ...
    "P-axis using mean amplitude: atan2(aVF, I)", ...
    "P-axis using mean amplitude: atan2(2*aVF, sqrt(3)*I)", ...
    "P-axis using mean amplitude: atan2(2*II-I, sqrt(3)*I)", ...
    ...
    "R-axis using mean amplitude: atan2(aVF, I)", ...
    "R-axis using mean amplitude: atan2(2*aVF, sqrt(3)*I)", ...
    "R-axis using mean amplitude: atan2(2*II-I, sqrt(3)*I)", ...
    ...
    "T-axis using median amplitude: atan2(aVF, I)", ...
    "T-axis using median amplitude: atan2(2*aVF, sqrt(3)*I)", ...
    "T-axis using median amplitude: atan2(2*II-I, sqrt(3)*I)", ...
    ...
    "P-axis using median amplitude: atan2(aVF, I)", ...
    "P-axis using median amplitude: atan2(2*aVF, sqrt(3)*I)", ...
    "P-axis using median amplitude: atan2(2*II-I, sqrt(3)*I)", ...
    ...
    "R-axis using median amplitude: atan2(aVF, I)", ...
    "R-axis using median amplitude: atan2(2*aVF, sqrt(3)*I)", ...
    "R-axis using median amplitude: atan2(2*II-I, sqrt(3)*I)" ...
    };    

end