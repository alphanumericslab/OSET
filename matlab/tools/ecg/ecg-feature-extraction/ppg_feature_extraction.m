function [ppg_features_vector, ppg_feature_info, ppg_fiducial_position, exit_flag] = ...
    ppg_feature_extraction(ppg_data, fs, ecg_rpeaks_index, lead_names, flatten_flag)

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
%   2. HRV features (7)
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


[C, T] = size(ppg_data); % number of channels

if nargin<3  || isempty(lead_names)
    for c = 1:C
        lead_names{c} = ['ch',num2str(c)];
    end
end


if nargin < 8 || isempty(flatten_flag)
    flatten_flag = true;
end

exit_flag = 0;


%% Feature Extraction


n_features = 12 ;

ppg_features_vector = [];
ppg_fiducial_position = cell(C,1);
if flatten_flag>0
    ecg_feature_names_ch = cell(1,C*n_features);
end

ppg_feature_names = [];

for c = 1:C

    data_channel = ppg_data(c, :);
    try

        flag_post_processing = 1;
        ppg_positions = fiducial_det_ppg(data_channel, ecg_rpeaks_index, fs, flag_post_processing);
        ppg_fiducial_position{c} = ppg_positions;

        ppg_sys_peaks_index = ppg_positions.onset;
        % Intervals
        data_rr_sys = em_interval_calc(ppg_sys_peaks_index'*1000/fs);
        if length(data_rr_sys)<length(ppg_sys_peaks_index)
            data_rr_sys(end+1) = data_rr_sys(end);
        end
        data_rr_sys(data_rr_sys>2000) = nan;
        data_rr_sys = fillmissing(data_rr_sys,'linear');


        % Morphological: single wave intervals
        dat_ppg.ppg_onset_to_systolic_peak = ([ppg_positions.sys_peak'-ppg_positions.onset']*1000/fs);
        dat_ppg.systole = ([ppg_positions.d_notch'-ppg_positions.onset']*1000/fs);

        temp_dias = ([ppg_positions.onset(2:end)'-ppg_positions.d_notch(1:end-1)']*1000/fs);
        temp_dias(end+1) = temp_dias(end);
        temp_dias(temp_dias<100) = 100;
        temp_dias(temp_dias>1000) = 1000;
        dat_ppg.diastole = temp_dias;

        obs_ppg{1,1} = dat_ppg.ppg_onset_to_systolic_peak';
        obs_ppg{2,1} = dat_ppg.systole';
        obs_ppg{3,1} = dat_ppg.diastole';
        obs_ppg{4,1} = data_rr_sys';

        for e = 1:size(obs_ppg,1)
            channels_observations{1,e} = obs_ppg{e,1}(2:end);
            channels_observations{2,e} = obs_ppg{e,1}(1:end-1);
        end

        for e = 1:size(obs_ppg,1)
            mn_ppg(1,e) =  mean(channels_observations{1,e},'omitnan');
            rmssd_ppg(1,e) =  sqrt(0.5)*std(channels_observations{1,e}-channels_observations{2,e},[],'omitnan');%/ std(channels_observations{1,e}+channels_observations{2,e});%mean(channels_observations{1,e}(this_block_ecg_min>0)) - mean(channels_observations{1,e}(this_block_ecg_min<=0));
            sdnn_ppg(1,e) =  sqrt(0.5)*std(channels_observations{1,e},[],'omitnan');
        end

        % Append to the overall feature vector
        all_features = [mn_ppg, rmssd_ppg, sdnn_ppg];
        if flatten_flag>0
            ppg_features_vector = cat(2, ppg_features_vector, all_features);
        else
            ppg_features_vector = cat(1, ppg_features_vector, all_features);
        end

    catch ME
        % Error handling: display message and set NaN values for problematic channel
        disp(ME.message);
        fprintf("Error in processing signal_channel: %s, %d\n", c);
        ppg_features_vector = cat(2, ppg_features_vector, nan(1, n_features));
        exit_flag = exit_flag -1;

    end

    if c==1 || isempty(ppg_feature_names)

        ppg_features_units = repmat({'ms'}, 1, size(all_features,2));

        % Define feature info
        ppg_feature_names = {'onsys_mean', 'systole_mean',  'diastole_mean', 'puls_mean',...
            'onsys_rmssd', 'systole_rmssd',  'diastole_rmssd', 'puls_rmssd', ...
            'onsys_sdnn', 'systole_sdnn',  'diastole_sdnn', 'puls_sdnn'};

        ppg_feature_description = {"Mean onset to sysolic peak intervals","Mean systolic interval","Mean diastolic interval","Mean pulse interval", ...
            "RMSSD of onset to sysolic peak intervals", "RMSSD of systolic interval","RMSSD of diastolic interval","RMSSD of pulse interval", ...
            "SDNN onset to sysolic peak intervals", "SDNN systolic interval","SDNN diastolic interval","SDNN pulse interval"};
    end

    if flatten_flag>0
        for n = 1:length(ppg_feature_names)
            ecg_feature_names_ch{(c-1)*n_features+n} = [ppg_feature_names{n},'_', lead_names{c}];
        end
    end

end

if flatten_flag==0
    ecg_feature_names_ch = ppg_feature_names;
else
    % Repeat cell array C times
    ppg_features_units = repmat(ppg_features_units, 1, C);
    ppg_feature_description = repmat(ppg_feature_description, 1, C);
end

ppg_feature_description = ppg_feature_description';
ecg_feature_names_ch = ecg_feature_names_ch';
ppg_features_units = ppg_features_units';

ppg_feature_info.names = ecg_feature_names_ch;
ppg_feature_info.units = ppg_features_units;
ppg_feature_info.description = ppg_feature_description;

end
