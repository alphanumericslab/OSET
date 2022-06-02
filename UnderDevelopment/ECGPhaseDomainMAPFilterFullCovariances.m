function [data_posterior_est, data_prior_est, n_var] = ECGPhaseDomainMAPFilterFullCovariances(data, peaks, params)
% An ECG denoiser based on a data driven MAP estimator in the phase domain
%
% Usage:
%   [data_posterior_est, data_prior_est, n_var] = ECGPhaseDomainMAPFilter(data, peaks, params)
%
% Inputs:
%   data: single or multichannel ECG signal with row-wise channels
%   peaks: R-peaks vector (common between all ECG channels)
%   params: a structure containing the algorithm parameters
%       params.bins: number of phase bins used for the ECG phase domain representation
%       params.BEAT_AVG_METHOD: Beat averaging method ('MEAN' or 'MEDIAN')
%       params.NOISE_VAR_EST_METHOD: noise variance estimation method
%           ('MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE')
%       params.avg_bins: Number of bins used in noise variance estimation
%           based on averaging ('AVGLOWER', 'MEDLOWER')
%       params.SMOOTH_PHASE: phase smoothing method ('BYPASS', 'MA', 'GAUSSIAN')
%       params.gaussianstd: smoothing window STD for params.SMOOTH_PHASE='GAUSSIAN'
%       params.wlen_time: smoothing time window for params.SMOOTH_PHASE='MA'
%       params.wlen_phase: smoothing phase window for params.SMOOTH_PHASE='MA'
%       params.nvar: stationary noise variance (if available). Estimated internally, if not given as input
%       params.nvar_factor: noise variance over/under estimation factor. Default = 1
%       params.plotresults: plot results or not (true or false)
%
% Outputs:
%   data_posterior_est: ECG signal based on posterior distribution
%       estimates (row-wise channels)
%   data_prior_est: ECG signal based on prior distribution estimates
%       (row-wise channels)
%   n_var: the noise variance estimate used in the algorithm (estimated internally)
%
% Refs:
%   1- ECG phase domain definition: R. Sameni, C. Jutten, and M. B. Shamsollahi.
%   "Multichannel electrocardiogram decomposition using periodic component
%   analysis." IEEE transactions on biomedical engineering 55.8 (2008): 1935-1940.
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

if ~isfield(params, 'BEAT_AVG_METHOD') % beat averaging method
    params.BEAT_AVG_METHOD = 'MEAN';
    warning(['Beat averaging method not defined; using ', params.BEAT_AVG_METHOD, ' as default']);
end

if ~isfield(params, 'nvar_factor') % noise variance over/under estimation factor
    params.nvar_factor = 1;
end


phaseshift = pi/params.bins;
pphase = PhaseShifting(phase, phaseshift); % phase shifting to compensate half a bin shift

[M, end_of_beat_indexes] = ECGPhaseToMatrix(pphase, params.bins); % Convert ECG phase into a (phase x time) matrix

data_prior_est = zeros(size(data));
data_posterior_est = zeros(size(data));
for ch = 1 : size(data, 1)
    x = data(ch, :); % the active channel

    % construct a stacked matrix representation of the ECG beats in the phased omain
    x_stacked = zeros(length(end_of_beat_indexes), params.bins);
    for ll = 1 : length(end_of_beat_indexes)
        if ll > 1
            beat_start = end_of_beat_indexes(ll-1) + 1;
        else
            beat_start = 1;
        end
        beat_end = end_of_beat_indexes(ll);
        x_stacked(ll, :) = x(beat_start : beat_end) * (diag(max(1, sum(M(:, beat_start : beat_end), 2))) \ M(:, beat_start : beat_end))';
        %         x_stacked(ll, :) = x(beat_start : beat_end) * M(:, beat_start : beat_end)';
    end
    ECG_mean = mean(x_stacked, 1);
    ECG_median = median(x_stacked, 1);

    [~, ~, meanphase, ~, ECGSamplesPerBin] = MeanECGExtraction(x, phase, params.bins, 1); % mean ECG extraction

    switch params.BEAT_AVG_METHOD
        case 'MEAN'
            ECG_avg = ECG_mean;
        case 'MEDIAN'
            ECG_avg = ECG_median;
    end

    %     x_stacked_zero_mean = x_stacked - ECG_avg(ones(1, size(x_stacked, 1)), :);
    %     K_x_ph = x_stacked_zero_mean' * x_stacked_zero_mean / size(x_stacked_zero_mean, 1);
    K_x_ph = cov(x_stacked); % identical to the above lines if called as: cov(x_stacked, 1)
    ECG_std = sqrt(diag(K_x_ph))';

    switch params.NOISE_VAR_EST_METHOD
        case 'MIN' % Method 1: min std
            noise_std_est = min(ECG_std);
        case 'AVGLOWER' % Method 2: average of the smallest std
            if ~isfield(params, 'avg_bins')
                params.avg_bins = 3;
                warning(['Undefined numer of averaging bins; using ', num2str(params.avg_bins) ' as default']);
            end
            [ECG_std_up_sorted, std_sorted_indexes] = sort(ECG_std);
            ECGSamplesPerBin_sorted = ECGSamplesPerBin(std_sorted_indexes);
            bns = 1 : params.avg_bins;
            noise_std_est = sqrt(sum(ECG_std_up_sorted(bns).^2 .* (ECGSamplesPerBin_sorted(bns) - 1)) / (sum(ECGSamplesPerBin_sorted(bns)) - 1)); % Note: recovers the original variances and renormalizes by total (N-1) to obtain an unbiased estimator
        case 'MEDLOWER' % Method 3: median of the smallest std
            if ~isfield(params, 'avg_bins')
                params.avg_bins = 3;
                warning(['Undefined numer of averaging bins; using ', num2str(params.avg_bins) ' as default']);
            end
            ECG_std_up_sorted = sort(ECG_std);
            noise_std_est = sqrt(median(ECG_std_up_sorted(1 : params.avg_bins).^2)); % Note: the sequence of sqrt->mean->square is to avoid swapping averaging and sqrt, which is nonlinear
        case 'PERCENTILE' % Method 4: percentiles
            if ~isfield(params, 'p')
                params.p = 0.5;
                warning(['Undefined averaging percentile; using ', num2str(params.p) ' as default']);
            end
            noise_std_est = prctile(ECG_std, params.p); % the desired percentile
        otherwise
            error('Undefined noise variance estimation method');
    end

    if ~isfield(params, 'nvar')
        n_var = params.nvar_factor * noise_std_est^2; % noise variance estimate
    else
        n_var = params.nvar;
    end

    % disp(['nvar estimate = ' num2str(sqrt(n_var))])

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
    if 0 % under dev.
        K_s_ph = K_x_ph - n_var * eye(params.bins);
        rho = K_s_ph / n_var;
        sig_len = length(x);
        T = M_smoothed / length(end_of_beat_indexes);
        %     K_s_tm = M' * K_s_ph * M;
        %     K_x_tm = M' * K_s_ph * M + n_var * eye(size(M, 2));
        data_prior_est(ch, :) = (0.5*length(end_of_beat_indexes)) * (T'*T) * x';%ECG_avg * T; % prior estimate based on average beat repetition over time
        % MAP estimate of each ECG sample assuming a Gaussian distribution for
        % the ECG samples and the noise (derived theoretically)
        %     data_posterior_est(ch, :) = ((T' * rho * T) * (eye(sig_len) - T' * pinv(pinv(rho) + (T * T')) * T) * (eye(sig_len) - (T'*T)) + (T'*T)) * x';
        data_posterior_est(ch, :) = (T' * rho * T) * (eye(sig_len) - T' * pinv(pinv(rho) + (T * T')) * T) * (x' - data_prior_est(ch, :)') + data_prior_est(ch, :)';
    else % under dev.
        K_s_ph = K_x_ph - n_var * eye(params.bins);

        % make matrix SPD
        [U,S,V] = svd(K_s_ph);
        S(S < 0) = 0;
        K_s_ph = U * S * V';

        rho = K_s_ph / n_var;
        sig_len = length(x);
        T = M_smoothed;
        %     K_s_tm = M' * K_s_ph * M;
        %     K_x_tm = M' * K_s_ph * M + n_var * eye(s ize(M, 2));
        %T2 = diag(sum(T, 2)) \ T;
        T2 = T / diag(sum(T, 1));
        
%         data_prior_est(ch, :) = T' * (T2 * x');%ECG_avg * T; % prior estimate based on average beat repetition over time
        data_prior_est(ch, :) = ECG_avg * T2;%ECG_avg * T; % prior estimate based on average beat repetition over time

        % MAP estimate of each ECG sample assuming a Gaussian distribution for
        % the ECG samples and the noise (derived theoretically)
        %     data_posterior_est(ch, :) = ((T' * rho * T) * (eye(sig_len) - T' * pinv(pinv(rho) + (T * T')) * T) * (eye(sig_len) - (T'*T)) + (T'*T)) * x';

        %         data_posterior_est(ch, :) = (T2' * rho * T2) * (eye(sig_len) - T2' * pinv(pinv(rho) + (T2 * T2')) * T2) * (x' - data_prior_est(ch, :)') + data_prior_est(ch, :)';

        Delta = x' - data_prior_est(ch, :)';
        term1 = T2 * Delta;
        term2 = T2' * rho;
        term3 = T2 * T2';
        term3_ = tril(triu(term3, -round(params.bins/2)), round(params.bins/2)); % removes the phase glithces at the corners of term3
        data_posterior_est(ch, :) = (term2 * term1) - term2 * term3_ * pinv(pinv(rho) + term3_) * term1 + data_prior_est(ch, :)';

        if 0
            figure
            plot(x)
            hold on
            plot(data_prior_est(ch, :));
            plot(data_posterior_est(ch, :));
            grid
        end
    end

    %     plot(x)
    %     hold on
    %     plot(data_prior_est(ch, :))
    %     plot(data_posterior_est(ch, :))
    %     grid

    if isfield(params, 'plotresults')
        if params.plotresults == true
            % if ch == 1
            %     figure
            %     plot(phase)
            %     hold on
            %     plot(pphase)
            %     grid
            %     legend('phase', 'pphase');
            % end

            figure
            errorbar(meanphase, ECG_avg, ECG_std/2);
            hold on
            plot(meanphase, ECG_avg);
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

