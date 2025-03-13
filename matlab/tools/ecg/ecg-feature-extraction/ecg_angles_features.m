function featureset = ecg_angles_features(data, position, fs)
    % featureset = ecg_angles_features(data, position, fs)
    % Extract features related to the ECG -Angles
    %
    % Inputs:
    %   data: ECG signal (1D array).
    %   position: fiducial points of ECG signal
    %   fs: Sampling frequency (Hz)
    %
    % Output:
    %   featureset: Stucture containing features related to the ECG -Angles
    %   (degree)
    %
    % Author:
    %   Seyedeh Somayyeh Mousavi 
    %   Emory University, Georgia, USA
    %   Email: bmemousavi@gmail.com
    %   Date: SEP 24, 2024
    % Author: Sajjad Karimi
    % Date: Mar 18, 2025

    %% Constant value
    convert_s_ms =1000;
    %% Define ECG points
    p = position.P;
    qrs_onset = position.QRSon;
    rpeak = position.R;
    qrs_offset = position.QRSoff;
    t = position.T;

    %% Calculate intervals
    rt_interval = convert_s_ms * ( t - rpeak) / fs;
    pr_interval = convert_s_ms * (rpeak - p) / fs;
    qr_interval = convert_s_ms * (rpeak - qrs_onset) / fs;
    rs_interval = convert_s_ms * (qrs_offset - rpeak) / fs;

    % Initialize rt_amp with NaNs of the same size as rpeak
    rt_amp = NaN * ones(size(rpeak)); 
    % Loop over each index and compute the difference
    for i = 1:length(rpeak)
          if ~isnan(rpeak(i)) && ~isnan(t(i)) 
                rt_amp(i) = data(rpeak(i)) - data(t(i));
          end
    end

    % Initialize pr_amp with NaNs of the same size as rpeak
    pr_amp = NaN * ones(size(rpeak)); 
    % Loop over each index and compute the difference
    for i = 1:length(rpeak)
          if ~isnan(rpeak(i)) && ~isnan(p(i)) 
                pr_amp(i) = data(rpeak(i)) - data(p(i));
          end
    end

    % Initialize qr_amp with NaNs of the same size as rpeak
    qr_amp = NaN * ones(size(rpeak)); 
    % Loop over each index and compute the difference
    for i = 1:length(rpeak)
          if ~isnan(rpeak(i)) && ~isnan(qrs_onset(i)) 
                qr_amp(i) = data(rpeak(i)) - data(qrs_onset(i));
          end
    end

    % Initialize rs_amp with NaNs of the same size as rpeak
    rs_amp = NaN * ones(size(rpeak)); 
    % Loop over each index and compute the difference
    for i = 1:length(rpeak)
          if ~isnan(rpeak(i)) && ~isnan(qrs_offset(i)) 
                rs_amp(i) = data(rpeak(i)) - data(qrs_offset(i));
          end
    end

    % Calculte tan
    ratio_rt = rt_amp./ rt_interval;
    ratio_pr = pr_amp./ pr_interval;
    ratio_qr = qr_amp./ qr_interval;
    ratio_rs = rs_amp./ rs_interval;

    % Omit NaN values
    valid_ratio_rt = ratio_rt(~isnan(ratio_rt));
    valid_ratio_pr = ratio_pr(~isnan(ratio_pr));
    valid_ratio_qr = ratio_qr(~isnan(ratio_qr));
    valid_ratio_rs = ratio_rs(~isnan(ratio_rs));

    %% Check if any valid values are present for valid_ratio_rt
    if isempty(valid_ratio_rt)
        mean_rt_angles = NaN;
        std_rt_angles = NaN;
        median_rt_angles = NaN;
    else

        valid_ratio_rt_angles = rad2deg(atan(valid_ratio_rt));
        mean_rt_angles = mean(valid_ratio_rt_angles);
        std_rt_angles = std(valid_ratio_rt_angles);
        median_rt_angles = median(valid_ratio_rt_angles);
    end
    %% Check if any valid values are present for valid_ratio_pr
    if isempty(valid_ratio_pr)
        mean_pr_angles = NaN;
        std_pr_angles = NaN;
        median_pr_angles = NaN;
    else

        valid_ratio_pr_angles = rad2deg(atan(valid_ratio_pr));
        mean_pr_angles = mean(valid_ratio_pr_angles);
        std_pr_angles = std(valid_ratio_pr_angles);
        median_pr_angles = median(valid_ratio_pr_angles);
    end
    %% Check if any valid values are present for valid_ratio_qr
    if isempty(valid_ratio_qr)
        mean_qr_angles = NaN;
        std_qr_angles = NaN;
        median_qr_angles = NaN;
    else

        valid_ratio_qr_angles = rad2deg(atan(valid_ratio_qr));
        mean_qr_angles = mean(valid_ratio_qr_angles);
        std_qr_angles = std(valid_ratio_qr_angles);
        median_qr_angles = median(valid_ratio_qr_angles);
    end

    %% Check if any valid values are present for valid_ratio_rs
    if isempty(valid_ratio_rs)
        mean_rs_angles = NaN;
        std_rs_angles = NaN;
        median_rs_angles = NaN;
    else

        valid_ratio_rs_angles = rad2deg(atan(valid_ratio_rs));
        mean_rs_angles = mean(valid_ratio_rs_angles);
        std_rs_angles = std(valid_ratio_rs_angles);
        median_rs_angles = median(valid_ratio_rs_angles);
    end

    %%  Results
    % Store pr_angles results
    featureset.mean_pr_angles = mean_pr_angles;
    featureset.std_pr_angles = std_pr_angles;
    featureset.median_pr_angles = median_pr_angles;

    % Store qr_angles results
    featureset.mean_qr_angles = mean_qr_angles;
    featureset.std_qr_angles = std_qr_angles;
    featureset.median_qr_angles = median_qr_angles;

    % Store rs_angles results
    featureset.mean_rs_angles = mean_rs_angles;
    featureset.std_rs_angles = std_rs_angles;
    featureset.median_rs_angles = median_rs_angles;

    % Store rt_angles results
    featureset.mean_rt_angles = mean_rt_angles;
    featureset.std_rt_angles = std_rt_angles;
    featureset.median_rt_angles = median_rt_angles;
end
