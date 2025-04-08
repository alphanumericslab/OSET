function [ppg_features_vector, ppg_feature_info, ppg_fiducial_position, exit_flag] = ...
    ppg_feature_extraction(ppg_data, fs, ecg_rpeaks_index, lead_names, flatten_flag)

% Description: Extract features from one record of multi-channel PPG signal
%
% INPUT:
% ppg_data - PPG signal as matrix CxT which C indicates number of PPG channels.
% fs   - Sampling frequency of the PPG signal (in Hz).
% ecg_rpeaks_index - Index of R peaks from ECG signal for PPG fiducial point detection
% lead_names - a cell array containing each PPG channel name
% flatten_flag - a boolian flag used to concatenate features across all channels,
%                either together (true) or underneath (false)  (default true)

% DEPENDENCIES:
% 1. OSET package

% OUTPUT:
% ppg_features_vector - includes the following features per each channel:
%   1. Mean onset to systolic peak intervals
%   2. Mean systolic intervals
%   3. Mean diastolic intervals
%   4. Mean pulse intervals
%   5. RMSSD of onset to systolic peak intervals
%   6. RMSSD of systolic intervals
%   7. RMSSD of diastolic intervals
%   8. RMSSD of pulse intervals
%   9. SDNN of onset to systolic peak intervals
%  10. SDNN of systolic intervals
%  11. SDNN of diastolic intervals
%  12. SDNN of pulse intervals
% ppg_feature_info - includes following information fields
%   ppg_feature_names - name of features for all channels
%   ppg_features_units - unit of all extracted features
%   ppg_fiducial_position - cell array of detected fiducial points of PPG
% exit_flag - with 0 indicates successful, and minus values (-K)
%   indicates error in processing K channels out of C channels

% Author:
%   Sajjad Karimi
%   Emory University, Georgia, USA
%   Email: sajjadkarimi91@gmail.com
%   Date: Mar 14, 2025


C = size(ppg_data,1); % number of channels

if nargin<4  || isempty(lead_names)
    for c = 1:C
        lead_names{c} = ['ch',num2str(c)];
    end
end


if nargin < 5 || isempty(flatten_flag)
    flatten_flag = 1;
end

exit_flag = 0;


%% Feature Extraction

n_features = 12;

% Initialize ppg_features_vector with NaN values based on flatten_flag
if flatten_flag > 0
    % For flattened features (all channels concatenated horizontally)
    ppg_features_vector = nan(1, C * n_features);
else
    % For stacked features (channels stacked vertically)
    ppg_features_vector = nan(C, n_features);
end

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
        T = length(ppg_sys_peaks_index);
        data_rr_sys = zeros(T,1);
        % Intervals
        data_rr_sys(1:T-1) = em_interval_calc(ppg_sys_peaks_index'*1000/fs);
        data_rr_sys(T) = data_rr_sys(T-1);

        data_rr_sys(data_rr_sys>2000) = nan;
        data_rr_sys = fillmissing(data_rr_sys,'linear');

        % Morphological: single wave intervals
        dat_ppg.ppg_onset_to_systolic_peak = ((ppg_positions.sys_peak'-ppg_positions.onset')*1000/fs);
        dat_ppg.systole = ((ppg_positions.d_notch'-ppg_positions.onset')*1000/fs);

        temp_dias = ppg_positions.onset';
        temp_dias(1:end-1) = ((ppg_positions.onset(2:end)'-ppg_positions.d_notch(1:end-1)')*1000/fs);
        temp_dias(end) = temp_dias(end-1);
        temp_dias(temp_dias<100) = 100;
        temp_dias(temp_dias>1000) = 1000;
        dat_ppg.diastole = temp_dias;

        obs_ppg = cell(4,1);
        obs_ppg{1,1} = dat_ppg.ppg_onset_to_systolic_peak';
        obs_ppg{2,1} = dat_ppg.systole';
        obs_ppg{3,1} = dat_ppg.diastole';
        obs_ppg{4,1} = data_rr_sys';

        mn_ppg = zeros(1,4);
        rmssd_ppg = zeros(1,4);
        sdnn_ppg = zeros(1,4);

        for e = 1:size(obs_ppg,1)
            mn_ppg(1,e) =  mean(obs_ppg{e,1},'omitnan');
            rmssd_ppg(1,e) =  sqrt(0.5)*std(diff(obs_ppg{e,1}),[],'omitnan');%/ std(channels_observations{1,e}+channels_observations{2,e});%mean(channels_observations{1,e}(this_block_ecg_min>0)) - mean(channels_observations{1,e}(this_block_ecg_min<=0));
            sdnn_ppg(1,e) =  sqrt(0.5)*std(obs_ppg{e,1},[],'omitnan');
        end

        % Append to the overall feature vector
        all_features = [mn_ppg, rmssd_ppg, sdnn_ppg];

        % Fill the pre-initialized ppg_features_vector based on flatten_flag
        if flatten_flag > 0
            % For flattened features (all channels concatenated horizontally)
            start_idx = (c-1) * n_features + 1;
            end_idx = c * n_features;
            ppg_features_vector(start_idx:end_idx) = all_features;
        else
            % For stacked features (channels stacked vertically)
            ppg_features_vector(c, :) = all_features;
        end

    catch ME
        % Error handling: display message and set NaN values for problematic channel
        disp(ME.message);
        fprintf("Error in processing signal_channel: %s, %d\n", c);
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
