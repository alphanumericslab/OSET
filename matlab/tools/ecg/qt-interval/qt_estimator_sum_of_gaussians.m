function fiducials = qt_estimator_sum_of_gaussians(ecg, fs, varargin)
% qt_estimator_sum_of_gaussians - Estimates Q-T interval and related parameters using ECG data.
%
% Syntax:
%   fiducials = qt_estimator_sum_of_gaussians(ecg, fs)
%   fiducials = qt_estimator_sum_of_gaussians(ecg, fs, params)
%
% Inputs:
%   ecg: Matrix of multichannel ECG data, where each row corresponds to a channel.
%   fs: Sampling frequency of the ECG data.
%   params (optional): A structure containing parameter settings:
%       w1: Window length of the first median filter used in baseline removal. Default: 0.75
%       w2: Window length of the second median filter used in baseline removal. Default: 0.9
%       f_hr: Heart rate factor for peak detection. Default: 1.2
%       amp_th: Amplitude threshold for peak detection. Default: 0.6
%       polarity_detection_mode: Polarity detection mode for peak detection. Default: 2
%       initial_q: Initial values for Q wave segments of interest. Default: [-0.080; -0.020]
%       initial_t: Initial values for T wave segments of interest. Default: [0.1; 0.5]
%       BayesianEstimation: perform Bayesian estimation or not (true/false).
%           Saves processing load if only ML estimation is required.
%           Default: true
%       snr: Signal-to-noise ratio for Bayesian approach. Default: 100
%       beta: Constant coefficient multiplied by the Q/T wave widths for
%           Q-T wave width determination. Default: 3. Can be a scalar such as
%           beta = 3 or a construct with beta.q = 2.1 and beta.t = 2.8
%       widest_qt: widest QT-interval. QT estimates beyond this threshold
%           are set to nan. Default: .6 seconds
%       earliest_q_onset: earliest Q-wave onset relative to the R-peak.
%           Q-wave onsets and offsets are both set to nan, if the Q onset
%           or offsets are estimated to be earlier. Default: 0.15 seconds
%       latest_t_offset: latest T-wave offset relative to the R-peak.
%           T-wave onsets and offsets are both set to nan, if the T onset
%           or offsets are estimated to be later. Default: 0.4 seconds
%       plot_results: Whether to plot results. Default: false
%
% Outputs:
%   fiducials: A structure containing various fiducial points and intervals
%       estimated using both ML and Bayesian approaches. The fiducials
%       structure contains the following fields, which are vectors of the
%       same length as the number of detected R-peaks:
%       - R_peaks: Detected R peak positions for each channel.
%       - QT_ML: QT intervals estimated using the Maximum Likelihood approach for each channel.
%       - Q_onsets_ML: Q wave onset positions estimated using the Maximum Likelihood approach for each channel.
%       - Q_offsets_ML: Q wave offset positions estimated using the Maximum Likelihood approach for each channel.
%       - T_onsets_ML: T wave onset positions estimated using the Maximum Likelihood approach for each channel.
%       - T_offsets_ML: T wave offset positions estimated using the Maximum Likelihood approach for each channel.
%       - QT_BYS: QT intervals estimated using the Bayesian approach for each channel.
%       - Q_onsets_BYS: Q wave onset positions estimated using the Bayesian approach for each channel.
%       - Q_offsets_BYS: Q wave offset positions estimated using the Bayesian approach for each channel.
%       - T_onsets_BYS: T wave onset positions estimated using the Bayesian approach for each channel.
%       - T_offsets_BYS: T wave offset positions estimated using the Bayesian approach for each channel.
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions on
%   Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
%   2023: First release
%
% Davood Fattahi and Reza Sameni 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Handle optional input arguments
if nargin > 2 && ~isempty(varargin{1})
    params = varargin{1};
else
    params = struct;
end

% Default parameter values
if isfield(params, 'w1')
    w1 = params.w1;
else
    w1 = 0.75; % window length of the first median filter used in baseline removal
end

if isfield(params, 'w2')
    w2 = params.w2;
else
    w2 = 0.9; % window length of the second median filter used in baseline removal
end

if isfield(params, 'f_hr')
    f_hr = params.f_hr;
else
    f_hr = 1.2;
end

if isfield(params, 'amp_th')
    amp_th = params.amp_th;
else
    amp_th = 0.6;
end

if isfield(params, 'polarity_detection_mode')
    polarity_detection_mode = params.polarity_detection_mode;
else
    polarity_detection_mode = 2;
end

if isfield(params, 'initial_q')
    soi_initial.q = params.initial_q;
else
    soi_initial.q = [-0.080; -0.020];
end

if isfield(params, 'initial_t')
    soi_initial.t = params.initial_t;
else
    soi_initial.t = [0.1; 0.5];
end

if isfield(params, 'BayesianEstimation')
    BayesianEstimation = params.BayesianEstimation;
else
    BayesianEstimation = true;
end

if isfield(params, 'snr')
    snr = params.snr;
else
    snr = 100;
end

if isfield(params, 'beta')
    beta = params.beta;
else
    beta = 3;
end

if ~isstruct(beta)
    Beta = beta;
    beta = struct;
    beta.q = Beta;
    beta.t = Beta;
end

if isfield(params, 'plot_results')
    plot_results = params.plot_results;
else
    plot_results = false;
end

if isfield(params, 'widest_qt')
    widest_qt = params.widest_qt;
else
    widest_qt = .6;
end

if isfield(params, 'earliest_q_onset')
    earliest_q_onset = params.earliest_q_onset;
else
    earliest_q_onset = .150;
end

if isfield(params, 'latest_t_offset')
    latest_t_offset = params.latest_t_offset;
else
    latest_t_offset = .400;
end


% Number of channels
num_ch = size(ecg, 1);

% Signal length
len = size(ecg, 2);

% Initialize fiducials structure
fiducials.R_peaks = cell(1, num_ch);
fiducials.QT_ML = cell(1, num_ch);
fiducials.Q_onsets_ML = cell(1, num_ch);
fiducials.Q_offsets_ML = cell(1, num_ch);
fiducials.T_onsets_ML = cell(1, num_ch);
fiducials.T_offsets_ML = cell(1, num_ch);
if BayesianEstimation
    fiducials.QT_BYS = cell(1, num_ch);
    fiducials.Q_onsets_BYS = cell(1, num_ch);
    fiducials.Q_offsets_BYS = cell(1, num_ch);
    fiducials.T_onsets_BYS = cell(1, num_ch);
    fiducials.T_offsets_BYS = cell(1, num_ch);
end

% Baseline removal using median filtering
ecg = ecg - baseline_sliding_window(baseline_sliding_window(ecg, round(w1 * fs), 'md'), round(w2 * fs), 'mn');

for ch = 1 : num_ch
    % Peak detection using amplitude threshold
    [~, r_peak_indexes] = peak_det_amp_threshold(ecg(ch, :), f_hr / fs, amp_th, polarity_detection_mode, 'MEDIAN');
    fiducials.R_peaks{ch} = r_peak_indexes;

    % Extract Q and T wave segments of interest
    [~, Q, ~, ~, T] = ecgWavesSoI_RR(ecg(ch, :), r_peak_indexes, 1);
    soi2.q = Q(:, [1, 3]) / fs;
    soi2.t = T(:, [1, 3]) / fs;

    [~, Q, ~, ~, T] = ecgWavesSoI_k(ecg(ch, :), r_peak_indexes, 1);
    soi3.q = Q(:, [1, 3]) / fs;
    soi3.t = T(:, [1, 3]) / fs;

    % Initial parameters and bounds
    p0 = [];
    lb = [];
    ub = [];

    % Optimization settings
    options = struct('SpecifyObjectiveGradient', true);

    % ML approach
    [mlGaussParams, ~, soi_ML, waveParams_ML, qtInt_ML] = qtParamsGausFit_cl(ecg(ch, :)', fs, r_peak_indexes, p0, soi_initial, lb, ub, beta, options);

    % Store ML results in the fiducials structure
    fiducials.QT_ML{ch} = qtInt_ML;
    fiducials.Q_onsets_ML{ch} = mlGaussParams.q(3, :) - beta.q * mlGaussParams.q(2, :);
    fiducials.Q_offsets_ML{ch} = mlGaussParams.q(3, :) + beta.q * mlGaussParams.q(2, :);
    fiducials.T_onsets_ML{ch} = mlGaussParams.t(3, :) - beta.t * mlGaussParams.t(2, :);
    fiducials.T_offsets_ML{ch} = mlGaussParams.t(3, :) + beta.t * mlGaussParams.t(2, :);

    % Bayesian approach
    varNoise = var(ecg(ch, :)) / 10.^(snr/10); % calculate noise variance from presumed or estimated SNR

    PrMu.q = mean(mlGaussParams.q(:, :), 2, "omitnan");
    PrCov.q = cov(mlGaussParams.q(:, :)', "omitrows");

    PrMu.t = mean(mlGaussParams.t(:, :), 2, "omitnan");
    PrCov.t = cov(mlGaussParams.t(:, :)', "omitrows");

    if BayesianEstimation
        [bysGaussParams, ~, soi_BYS, waveParams_BYS, qtInt_BYS] = qtParamsGausFit_cl(ecg(ch, :)', fs, r_peak_indexes, p0, soi_ML, lb, ub, beta, PrMu, PrCov, varNoise);

        % Store Bayesian results in the fiducials structure
        fiducials.QT_BYS{ch} = qtInt_BYS;
        fiducials.Q_onsets_BYS{ch} = bysGaussParams.q(3, :) - beta.q * bysGaussParams.q(2, :);
        fiducials.Q_offsets_BYS{ch} = bysGaussParams.q(3, :) + beta.q * bysGaussParams.q(2, :);
        fiducials.T_onsets_BYS{ch} = bysGaussParams.t(3, :) - beta.t * bysGaussParams.t(2, :);
        fiducials.T_offsets_BYS{ch} = bysGaussParams.t(3, :) + beta.t * bysGaussParams.t(2, :);
    end

    if plot_results
        time = (1 : size(ecg, 2))/fs; % time stamp

        % Plot the evaluated gaussians on the signal
        figure; p1 = plot(time, ecg(ch, :)); % plot the channel
        hold on

        for j = 1 : length(r_peak_indexes)
            tq = soi_ML.q(j, 1) : 1/fs : soi_ML.q(j, 2);
            tt = soi_ML.t(j, 1) : 1/fs : soi_ML.t(j, 2);
            p2 = plot(tq + r_peak_indexes(j)/fs, GausVal(tq, mlGaussParams.q(:, j)), 'r-');
            p3 = plot(tt + r_peak_indexes(j)/fs, GausVal(tt, mlGaussParams.t(:, j)), 'r-');
            q_onset_index = round(fs * waveParams_ML.q(1, j) + r_peak_indexes(j));
            t_offset_index = round(fs * waveParams_ML.t(2, j) + r_peak_indexes(j));
            if q_onset_index > 0 && q_onset_index <= len && t_offset_index > 0 && t_offset_index <= len
                p4 = plot(time([q_onset_index, t_offset_index]), ecg(ch, [q_onset_index, t_offset_index]), 'c*');
            end
        end
        if p4
            legend([p1, p2, p3, p4], 'ecg', 'Q-segment model', 'T-segment model', 'q/t onset/offset')
        else
            legend([p1, p2, p3], 'ecg', 'Q-segment model', 'T-segment model')
        end
        clear p4
        title 'ML framework'

        if BayesianEstimation
            % Plot the evaluated gaussians on the signal
            figure; p1 = plot(time, ecg(ch, :)); % plot the channel
            hold on

            for j = 1 : length(r_peak_indexes)
                tq = soi_BYS.q(j, 1) : 1 / fs : soi_BYS.q(j, 2);
                tt = soi_BYS.t(j, 1) : 1 / fs : soi_BYS.t(j, 2);
                p2 = plot(tq + r_peak_indexes(j) / fs, GausVal(tq, bysGaussParams.q(:, j)), 'r-');
                p3 = plot(tt + r_peak_indexes(j) / fs, GausVal(tt, bysGaussParams.t(:, j)), 'r-');
                q_onset_index = round(fs * waveParams_BYS.q(1, j) + r_peak_indexes(j));
                t_offset_index = round(fs * waveParams_BYS.t(2, j) + r_peak_indexes(j));
                if q_onset_index > 0 && q_onset_index <= len && t_offset_index > 0 && t_offset_index <= len
                    p4 = plot(time([q_onset_index, t_offset_index]), ecg(ch, [q_onset_index, t_offset_index]), 'c*');
                end
            end
            if p4
                legend([p1, p2, p3, p4], 'ecg', 'Q-segment model', 'T-segment model', 'q/t onset/offset')
            else
                legend([p1, p2, p3], 'ecg', 'Q-segment model', 'T-segment model')
            end
            title 'BYS framework'
        end
    end
end

fiducials = clean_estimates(fiducials, widest_qt, earliest_q_onset, latest_t_offset); % replace outliers with nans
end

function fiducials = clean_estimates(fiducials, widest_qt, earliest_q_onset, latest_t_offset)
% A function to replace out of range parameters with nans
for ch = 1 : length(fiducials)
    % Clean the QT interval
    if isfield(fiducials, 'QT_ML')
        for k = 1 : length(fiducials.QT_ML{ch})
            if fiducials.QT_ML{ch}(k) > widest_qt
                fiducials.QT_ML{ch}(k) = nan;
            end
        end
    end

    % Clean the Q onsets/offsets (ML approach)
    if isfield(fiducials, 'Q_onsets_ML') && isfield(fiducials, 'Q_offsets_ML')
        for k = 1 : length(fiducials.Q_onsets_ML{ch})
            if fiducials.Q_onsets_ML{ch}(k) < -earliest_q_onset || fiducials.Q_offsets_ML{ch}(k) < -earliest_q_onset
                fiducials.Q_onsets_ML{ch}(k) = nan;
                fiducials.Q_offsets_ML{ch}(k) = nan;
            end
        end
    end

    % Clean the T offsets/offsets (ML approach)
    if isfield(fiducials, 'T_onsets_ML') && isfield(fiducials, 'T_offsets_ML')
        for k = 1 : length(fiducials.T_offsets_ML{ch})
            if fiducials.T_onsets_ML{ch}(k) > latest_t_offset || fiducials.T_offsets_ML{ch}(k) > latest_t_offset
                fiducials.T_onsets_ML{ch}(k) = nan;
                fiducials.T_offsets_ML{ch}(k) = nan;
            end
        end
    end

    % Clean the Q onsets/offsets (Bayesian approach)
    if isfield(fiducials, 'Q_onsets_BYS') && isfield(fiducials, 'Q_offsets_BYS')
        for k = 1 : length(fiducials.Q_onsets_BYS{ch})
            if fiducials.Q_onsets_BYS{ch}(k) < -earliest_q_onset || fiducials.Q_offsets_BYS{ch}(k) < -earliest_q_onset
                fiducials.Q_onsets_BYS{ch}(k) = nan;
                fiducials.Q_offsets_BYS{ch}(k) = nan;
            end
        end
    end

    % Clean the T offsets/offsets (Bayesian approach)
    if isfield(fiducials, 'T_onsets_BYS') && isfield(fiducials, 'T_offsets_BYS')
        for k = 1 : length(fiducials.T_offsets_BYS{ch})
            if fiducials.T_onsets_BYS{ch}(k) > latest_t_offset || fiducials.T_offsets_BYS{ch}(k) > latest_t_offset
                fiducials.T_onsets_BYS{ch}(k) = nan;
                fiducials.T_offsets_BYS{ch}(k) = nan;
            end
        end
    end
end
end
