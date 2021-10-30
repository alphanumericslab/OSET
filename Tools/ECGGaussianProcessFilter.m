function [data_posterior_est, data_prior_est] = ECGGaussianProcessFilter(data, peaks, params)
% Description: An ECG denoiser based on a data driven Gaussian Process model
%   of the ECG samples and noise
%
% Usage:
%   [data_posterior_est, data_prior_est] = ECGGaussianProcessFilter(data, peaks, params)
%
% Inputs:
%   data: single or multichannel ECG signal with row-wise channels
%   peaks: R-peaks vector (common between all ECG channels)
%   params: a structure containing the algorithm parameters
%       params.bins: number of phase bins used for the ECG phase domain representation
%       params.NOISE_VAR_EST_METHOD: noise variance estimation method
%           ('MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE')
%       params.avg_bins: Number of bins used in noise variance estimation
%           based on averaging ('AVGLOWER', 'MEDLOWER')
%       params.SMOOTH_PHASE: phase smoothing method ('BYPASS', 'MA', 'GAUSSIAN')
%       params.gaussianstd: smoothing window STD for params.SMOOTH_PHASE='GAUSSIAN'
%       params.wlen_time: smoothing time window for params.SMOOTH_PHASE='MA'
%       params.wlen_phase: smoothing phase window for params.SMOOTH_PHASE='MA'
%       params.plotresults: plot results or not (true or false)
%
% Outputs:
%   data_posterior_est: ECG signal based on posterior distribution
%       estimates (row-wise channels)
%   data_prior_est: ECG signal based on prior distribution estimates
%       (row-wise channels)
%
% Copyright Reza Sameni, Oct 2021
% The Open-Source Electrophysiological Toolbox
% (https://github.com/alphanumericslab/OSET)

[phase, ~] = PhaseCalculation(peaks); % phase calculation

if ~isfield(params, 'bins') % number of phase bins
    params.bins = median(diff(find(peaks))); % use median of RR-intervals in samples by default
end

if ~isfield(params, 'NOISE_VAR_EST_METHOD') % noise variance estimation method
    params.NOISE_VAR_EST_METHOD = 'MIN'; % Min variance across average beat, by default
end

if ~isfield(params, 'SMOOTH_PHASE') % phase smoothing method
    params.SMOOTH_PHASE = 'BYPASS'; % no smoothing, by default
end

if ~isfield(params, 'n_mean') % noise average
    params.n_mean = 0; % zero, by default
end

phaseshift = pi/params.bins;
pphase = PhaseShifting(phase, phaseshift); % phase shifting to compensate half a bin shift

% figure
% plot(phase)
% hold on
% plot(pphase)
% grid
% legend('phase', 'pphase');

M = ECGPhaseToMatrix(pphase, params.bins); % Convert ECG phase into a (phase x time) matrix
data_prior_est = zeros(size(data));
data_posterior_est = zeros(size(data));
for ch = 1 : size(data, 1)
    x = data(ch, :); % the active channel
    
    [ECG_mean, ECG_std, meanphase] = MeanECGExtraction(x, phase, params.bins, 1); % mean ECG extraction
    
    switch params.NOISE_VAR_EST_METHOD
        case 'MIN' % Method 1: min std
            noise_std_est = min(ECG_std);
        case 'AVGLOWER' % Method 2: average of the smallest std
            if ~isfield(params, 'avg_bins')
                params.avg_bins = 3;
                warning(['Undefined numer of averaging bins; using ', num2str(params.avg_bins) ' as default']);
            end
            ECG_std_up_sorted = sort(ECG_std);
            noise_std_est = sqrt(mean(ECG_std_up_sorted(1 : params.avg_bins).^2));
        case 'MEDLOWER' % Method 3: median of the smallest std
            if ~isfield(params, 'avg_bins')
                params.avg_bins = 3;
                warning(['Undefined numer of averaging bins; using ', num2str(params.avg_bins) ' as default']);
            end
            ECG_std_up_sorted = sort(ECG_std);
            noise_std_est = sqrt(median(ECG_std_up_sorted(1 : params.avg_bins).^2));
        case 'PERCENTILE' % Method 4: percentiles
            if ~isfield(params, 'p')
                params.p = 0.5;
                warning(['Undefined averaging percentile; using ', num2str(params.p) ' as default']);
            end
            noise_std_est = prctile(ECG_std, params.p); % the desired percentile
        otherwise
            error('Undefined noise variance estimation method');
    end
    n_var = noise_std_est^2; % noise variance estimate
    
    ECG_intrinsic_var = max(0, ECG_std.^2 - n_var); % average beat variance estimate
    
    % smooth the phase-time mapping matrix in time and phase
    switch params.SMOOTH_PHASE
        case 'BYPASS'
            M_smoothed = M;
        case 'MA'
            if ~isfield(params, 'wlen_time')
                params.wlen_time = 3;
                warning(['Undefined time smoothing window length; using ', num2str(params.wlen_time) ' as default']);
            end
            if ~isfield(params, 'wlen_phase')
                params.wlen_phase = 3;
                warning(['Undefined phase smoothing window length; using ', num2str(params.wlen_phase) ' as default']);
            end
            M_smoothed = filter2(ones(params.wlen_phase, params.wlen_time)/(params.wlen_phase * params.wlen_time), M);
        case 'GAUSSIAN'
            if ~isfield(params, 'gaussianstd')
                params.gaussianstd = 3;
                warning(['Undefined gaussian smoothing window length; using ', num2str(params.gaussianstd) ' as default']);
            end
            M_smoothed = imgaussfilt(M, params.gaussianstd);
    end
    
    % reconstruct the ECG
    data_prior_est(ch, :) = ECG_mean * M_smoothed; % prior estimate based on average beat repetition over time
    s_var = ECG_intrinsic_var * M_smoothed; % ECG beat variance estimate repeated over time
    
    % MAP estimate of each ECG sample assuming a Gaussian distribution for
    % the ECG samples and the noise (derived theoretically)
    data_posterior_est(ch, :) = (s_var .* (x - params.n_mean) + n_var * data_prior_est(ch, :)) ./ (s_var + n_var);
    
    if isfield(params, 'plotresults')
        if params.plotresults == true
            figure
            errorbar(meanphase, ECG_mean, ECG_std/2);
            hold on
            plot(meanphase, ECG_mean);
            grid
            legend('Errorbar', 'Average ECG');
            title('Average ECG beat \pm STD');
            xlabel('Phase (rad)');
            ylabel('Amplitude');
            set(gca, 'fontsize', 16)
            
            figure
            plot(meanphase, ECG_std);
            hold on
            plot(meanphase, sqrt(ECG_intrinsic_var));
            plot(meanphase, sqrt(n_var)*ones(1, params.bins));
            grid
            legend('Signal+Noise STD', 'Signal STD', 'Noise STD');
            title('Signal and noise STDs');
            xlabel('Phase (rad)');
            ylabel('Amplitude');
            set(gca, 'fontsize', 16)
            
            figure
            plot(sort(ECG_std));
            grid
            title('Sample STDs across the average beat (sorted)');
            xlabel('Sorted bin index');
            ylabel('STDd');
            set(gca, 'fontsize', 16)
            
            figure
            plot(x)
            hold on
            plot(x, 'linewidth', 2)
            plot(data_prior_est(ch, :), 'linewidth', 2)
            plot(data_posterior_est(ch, :), 'linewidth', 2)
            grid
            legend('Noisy ECG', 'Prior ECG estimate', 'Posterior ECG estimate')
            xlabel('time (samples)');
            ylabel('Amplitude');
            title('Filtering results');
            set(gca, 'fontsize', 16)
        end
    end
end

