% SAMPLE_DELI  Demo of LSIM-based ECG delineation with ecg_delineate_lsim
%
%   Runs ecg_delineate_lsim on every .mat record in this folder, plots
%   the detected fiducial points against the expert annotations and prints
%   the mean absolute error (ms) per fiducial point.
%
%   Each .mat record is expected to contain:
%       ecg           - ECG samples (T x channels); the first channel is used
%       fs            - sampling rate (Hz)
%       t_second      - time axis (s)
%       true_position - struct of expert annotations in samples
%                       (R, QRSon, QRSoff, Ton, T, Toff; NaN if absent)
%
%   Sajjad Karimi, Reza Sameni  2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

clear;
close all;
clc;

% update with path to .mat data files
db_folder = 'path/to/data';

local_db_files = dir(fullfile(db_folder, '*.mat')); % list of all mat files

% LSIM-Deli parameters (empty -> default value)
flag_post_processing = 1;           % use RR-interval priors to refine LSIM fiducials
flag_prune_P         = 0;           % prune low-quality P-waves using P_score
win_qrs              = [];          % feature window (s) for QRS on/off, default 0.01
win_T                = [];          % feature window (s) for T on/off, default 0.02
win_P                = [];          % feature window (s) for P on/off, default 0.02
max_clusters         = [];          % max FCM clusters, default adaptive (3..6)
twave_shape          = 'none';      % 'none', 'bi-phasic', 'min', 'max'
time_prior_mode      = 'disabled';  % 'disabled' or 'normalized'

fiducial_names = {'QRSon', 'QRSoff', 'Ton', 'T', 'Toff'};

%%

for m = 1:length(local_db_files)

    tic
    % load data and expert annotations
    in_fname = local_db_files(m).name(1:end-4);
    fprintf('[%d/%d] %s\n', m, length(local_db_files), in_fname);
    load(fullfile(db_folder, [in_fname, '.mat']), 'ecg', 'fs', 't_second', 'true_position');

    % ECG preprocessing
    ecg_denoised = ecg(:,1)';

    % NOTCH FILTERING THE ECG
    fc = 50.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % notch filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_denoised = filtfilt(b, a, ecg_denoised); % zero-phase non-causal filtering
    ecg_denoised = ecg_denoised - movmean(movmedian(ecg_denoised,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);
    ecg_denoised = lp_filter_zero_phase(ecg_denoised, 30/fs);


    % R-peak detection
    peak_detector_params.RETURN_SIGNAL_PEAKS = true; % return signal peaks or energy envelope peaks
    peak_detector_params.PLOT_RESULTS = false; % plot the results using the internal plot function of peak_det_likelihood or not
    peak_detector_params.PLOT_DIAGNOSTIC = false; % diagnostic mode (do not activate unless diving deep into the code! run only on short segments, since many figures are created)
    peak_detector_params.verbose = false; % reports all the default values for the internal parameters of peak_det_likelihood, which can be modified through this data structure if needed.
    peak_detector_params.REFINE_PEAKS = true;
    overlap_time = 1.0; % overlap between segments for continuity (1.0-2.0 seconds is enough)
    seg_len_time = 10.0; % segment length in seconds

    [~, ~, ecg_rpeaks_index] = peak_det_likelihood_long_recs(ecg_denoised, fs, seg_len_time, overlap_time, peak_detector_params);
    % To evaluate delineation independently of R-peak detection, use the
    % expert R-peaks instead:  ecg_rpeaks_index = true_position.R;
    % Passing [] lets ecg_delineate_lsim detect the R-peaks internally.

    % LSIM-Deli
    [lsim_positions, EXITFLAG] = ecg_delineate_lsim(ecg_denoised, fs, ecg_rpeaks_index, ...
        flag_post_processing, flag_prune_P, win_qrs, win_T, win_P, ...
        max_clusters, twave_shape, time_prior_mode);
    if strcmp(EXITFLAG.status, 'failed')
        warning('LSIM-Deli failed on %s: %s', in_fname, EXITFLAG.message.message);
    end

    % Mean absolute error (ms) against expert annotations, matched to the
    % nearest detected point within 150 ms
    for k = 1:length(fiducial_names)
        fn = fiducial_names{k};
        if ~isfield(true_position, fn), continue; end
        ref = true_position.(fn)(:);  ref(isnan(ref)) = [];
        det = lsim_positions.(fn)(:); det(isnan(det)) = [];
        if isempty(ref) || isempty(det), continue; end
        err_ms = 1000*min(abs(ref - det'), [], 2)/fs;
        err_ms(err_ms > 150) = [];
        fprintf('   %-7s MAE = %6.2f ms  (%d / %d matched)\n', fn, mean(err_ms), length(err_ms), length(ref));
    end
    toc

    % plot the results
    figure('Position', [130 130 1500 800]);
    lg = {};
    plot(t_second, ecg_denoised, LineWidth=1.5); lg = cat(2, lg, {'ECG'});
    hold on

    % {field, marker, label}
    lsim_markers = {'R','b*','R'; 'Pon','k+','Pon'; 'P','k*','P'; 'Poff','kx','Poff'; ...
        'QRSon','rx','QRSon'; 'QRSoff','rx','QRSoff'; 'Ton','m+','Ton'; 'T','m*','T'; 'Toff','mx','Toff'};
    expert_markers = {'QRSon','go','EXPERT-QRSon'; 'QRSoff','go','EXPERT-QRSoff'; ...
        'Ton','g+','EXPERT-Ton'; 'T','g*','EXPERT-T'; 'Toff','gx','EXPERT-Toff'};

    for k = 1:size(lsim_markers, 1)
        idx = lsim_positions.(lsim_markers{k,1}); idx(isnan(idx)) = [];
        if isempty(idx), continue; end
        plot(t_second(idx), ecg_denoised(idx), lsim_markers{k,2}, MarkerSize=12, LineWidth=2); lg = cat(2, lg, lsim_markers(k,3));
    end
    for k = 1:size(expert_markers, 1)
        if ~isfield(true_position, expert_markers{k,1}), continue; end
        idx = true_position.(expert_markers{k,1}); idx(isnan(idx)) = [];
        if isempty(idx), continue; end
        plot(t_second(idx), ecg_denoised(idx), expert_markers{k,2}, MarkerSize=10, LineWidth=2); lg = cat(2, lg, expert_markers(k,3));
    end

    grid on
    legend(lg, 'Interpreter', 'latex', 'orientation', 'horizontal', 'FontSize', 14)
    xlabel('time (sec)', Interpreter='latex', FontSize=14)
    title(in_fname, 'Interpreter', 'none')

end
