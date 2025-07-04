function [ecg_features_vector, ecg_feature_info, ecg_fiducial_position, exit_flag] = ...
    ecg_feature_extraction(ecg_data, fs, lead_names, n_svd, num_morpho_samples, norm_flag, feature_list, flatten_flag, sync_rpeaks)

% Description: Extract features from one record of multi-channel ECG signal
%
% INPUT:
% ecg_data - ECG signal as matrix CxT which C indicates number of ECG channels (in microvolts).
% fs   - Sampling frequency of the ECG signal (in Hz).
% lead_names - a cell array containing each ECG lead name
% n_svd - Number of eigenvalues for SVD features (default 5).
% num_morpho_samples - Number of samples for ECG average beat (default 25).
% norm_flag - A boolian value for normalizing ECG average beat before sampling (default true).
% feature_list - A cell array contaning a subset of {'snr', 'hrv','interval', 'amplitude', 'angle', 'svd', 'morpho'}
%                to return selected subset of all ECG features (default all features)
% flatten_flag - a boolian flag used to concatenate features across all channels,
%                either together (true) or underneath (false)  (default true)


% DEPENDENCIES:
% 1. OSET package

%
% OUTPUT:
% ecg_features_vector - includes the following features per each channel:
%   1. SNR features (2)
%   2. HRV features (11)
%   3. Features related to time interval (29)
%   4. Features related to amplitude and area under curve (32)
%   5. Features related to angle (12)
%   6. Complexity analysis features (3)
%   7. SVD features (default 5)
%   7. Average morphology (default 25)
% vecg_feature_info - includes following information feilds
%   ecg_feature_names - name of features for all channels
%   ecg_features_units - unit of all extracted features
%   ecg_fiducial_position - cell array of detected fiducial points of ECG
% exit_flag - with 0 indicates successful, and minus values (-K)
%   indicates error in processing K channels out of C channels

% Author: Seyedeh Somayyeh Mousavi
% Date: Jan 11, 2025
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com
% Author:
%   Sajjad Karimi
%   Emory University, Georgia, USA
%   Email: sajjadkarimi91@gmail.com
%   Date: Mar 14, 2025


[C, T] = size(ecg_data); % number of channels

if nargin<3  || isempty(lead_names)
    for c = 1:C
        lead_names{c} = ['ch',num2str(c)];
    end
end

if nargin<4 || isempty(n_svd)
    n_svd = 5;
end

if nargin < 5 || isempty(num_morpho_samples)
    num_morpho_samples = 25;
end

if nargin < 6 || isempty(norm_flag)
    norm_flag = true;
end


if nargin < 7 || isempty(feature_list)
    feature_list = {'snr', 'hrv','interval', 'amplitude', 'angle', 'hjorth', 'svd', 'morpho'};
end

if nargin < 8 || isempty(flatten_flag)
    flatten_flag = true;
end

if nargin < 9 || isempty(sync_rpeaks)
    sync_rpeaks = true;
end

exit_flag = 0;


%% Feature Extraction

feature_list_all = {'snr', 'hrv','interval', 'amplitude', 'angle', 'hjorth', 'svd', 'morpho'};
[~ , ilocb] = ismember(feature_list,feature_list_all);

num_features =         [2     11     29          32          12        3      n_svd   num_morpho_samples]; % per channel

n_features = sum(num_features(ilocb));

ecg_features_vector = [];
ecg_fiducial_position = cell(C,1);
if flatten_flag>0
    ecg_feature_names_ch = cell(1,C*n_features);
end

seg_len_time = min(10,T/fs); % 10 seconds or signal length for data with shorter than 10s
pad_len_time = 1;

ecg_feature_names = [];

if sync_rpeaks>0
    params.RETURN_SIGNAL_PEAKS = false;
    % Detect peaks based on recording length
    ecg_data = ecg_data./std(ecg_data,[],2);
    [~, rpeak_indexes_in] = peak_det_likelihood_long_recs(ecg_data, fs, seg_len_time, pad_len_time, params);
end

for c = 1:C

    data_channel = ecg_data(c, :);
    try


        if sync_rpeaks==0
            % Detect peaks based on recording length
            [~, rpeak_indexes_in] = peak_det_likelihood_long_recs(data_channel, fs, seg_len_time, pad_len_time);
        end
        % Run ECG fiducial points detector
        % position = fiducial_det_lsim(data_channel, R_peaks_indexes, fs);
        % Run wavedet_3D
        % heasig.nsig = 1;
        % heasig.freq = fs;
        % heasig.nsamp = length(data_channel);
        % [position, ~, ~] = wavedet_3D(data_channel', R_peaks_indexes, heasig);

        flag_post_processing = 1;
        position = fiducial_det_lsim( data_channel', rpeak_indexes_in, fs, flag_post_processing);

        ecg_fiducial_position{c} = position;
        rpeak_indexes = position.R;

        [featureset{1}, this_feature_info{1}, mean_beat] = ecg_snr_features(data_channel, rpeak_indexes);
        [featureset{2}, this_feature_info{2}] = ecg_hrv_features(rpeak_indexes, fs);
        [featureset{3}, this_feature_info{3}] = ecg_time_intervals_features(rpeak_indexes, position, fs);
        [featureset{4}, this_feature_info{4}] = ecg_area_amp_features(data_channel, position, fs);
        [featureset{5}, this_feature_info{5}] = ecg_angles_features(data_channel, position, fs);
        [featureset{6}, this_feature_info{6}] = hjorth_time_features(mean_beat,fs);
        [featureset{7}, this_feature_info{7}] = ecg_svd_features(data_channel, rpeak_indexes, n_svd);
        [featureset{8}, this_feature_info{8}] = morphology_features(mean_beat, num_morpho_samples, norm_flag);


        % Append to the overall feature vector
        all_features = cell2mat(featureset(ilocb));
        if flatten_flag>0
            ecg_features_vector = cat(2, ecg_features_vector, all_features);
        else
            ecg_features_vector = cat(1, ecg_features_vector, all_features);
        end
        this_feature_info = this_feature_info(ilocb);

        error_flag = 0;
    catch ME
        % Error handling: display message and set NaN values for problematic channel
        disp(ME.message);
        fprintf("Error in processing signal_channel: %s, %d\n", c);
        if flatten_flag>0
            ecg_features_vector = cat(2, ecg_features_vector, nan(1, n_features));
        else
            ecg_features_vector = cat(1, ecg_features_vector, nan(1, n_features));
        end

        exit_flag = exit_flag -1;
        error_flag = 1;

    end

    if (c==1 || isempty(ecg_feature_names)) && error_flag==0
        ecg_features_units = [];
        ecg_feature_description = [];

        for f = 1:length(this_feature_info)
            ecg_feature_names = [ecg_feature_names, this_feature_info{f}.names];
            ecg_features_units = [ecg_features_units, this_feature_info{f}.units];
            ecg_feature_description = [ecg_feature_description, this_feature_info{f}.description];
        end

    end
end

if exit_flag == -C
    ecg_feature_info = [];
    return
end

for c = 1:C
    if flatten_flag>0
        for n = 1:length(ecg_feature_names)
            ecg_feature_names_ch{(c-1)*n_features+n} = [ecg_feature_names{n},'_', lead_names{c}];
        end
    end
end

if flatten_flag==0
    ecg_feature_names_ch = ecg_feature_names;
else
    % Repeat cell array C times
    ecg_features_units = repmat(ecg_features_units, 1, C);
    ecg_feature_description = repmat(ecg_feature_description, 1, C);

end

ecg_feature_description = ecg_feature_description';
ecg_feature_names_ch = ecg_feature_names_ch';
ecg_features_units = ecg_features_units';

ecg_feature_info.names = ecg_feature_names_ch;
ecg_feature_info.units = ecg_features_units;
ecg_feature_info.description = ecg_feature_description;

end
