function [data_posterior_est, data_prior_est, n_var] = ecg_den_time_domain_gp(data, peaks, params)
% ecg_den_time_domain_gp - An ECG denoiser based on a data-driven Gaussian Process MAP estimator in the time domain
%
% Usage:
%   [data_posterior_est, data_prior_est, n_var] = ecg_den_time_domain_gp(data, peaks, params)
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
%       params.avg_beat_shaping_window: window length used to shape the average beat with a hamming window.
%       params.plotresults: plot results or not (true or false)
%
% Outputs:
%   data_posterior_est: ECG signal based on posterior distribution
%       estimates (row-wise channels)
%   data_prior_est: ECG signal based on prior distribution estimates
%       (row-wise channels)
%   n_var: the noise variance estimate used in the algorithm (estimated internally)
%
%   Revision History:
%       2021: First release
%       2023: Renamed from deprecated version ECGTimeDomainMAPFilter()
%
%   Reza Sameni, 2021-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% if ~isfield(params, 'bins') % number of phase bins
%     params.bins = median(diff(find(peaks))); % use median of RR-intervals in samples by default
% end

if ~isfield(params, 'NOISE_VAR_EST_METHOD') % noise variance estimation method
    params.NOISE_VAR_EST_METHOD = 'MIN'; % Min variance across average beat, by default
end

% if ~isfield(params, 'SMOOTH_PHASE') % phase smoothing method
%     params.SMOOTH_PHASE = 'BYPASS'; % no smoothing, by default
% end

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

SignalLen = size(data, 2);
peak_indexes = find(peaks);
NumBeats = length(peak_indexes);
event_width = max(diff([1, peak_indexes, SignalLen])) + 1; % Add the first and last indexes
if(mod(event_width, 2) == 0)
    event_width = event_width + 1;
end
half_len = floor(event_width/2); % Half the window length

% Make the time to average beat mapping matrix
M = zeros(event_width, SignalLen + 2 * event_width);
for kk = 1 : NumBeats
    cols = peak_indexes(kk) + (-half_len : half_len) + event_width;
    M(:, cols) = M(:, cols) + eye(event_width);
end
M = M(:, event_width + 1 : event_width + SignalLen);

% fill the beginning and end points (before the first and after the last beat)
first_set_index = peak_indexes(1) - half_len;
for m = 1 : first_set_index - 1
    M(:, m) = M(end:-1:1, 2*first_set_index - m - 1); % Repeat the pattern after the first beat
end
last_set_index = peak_indexes(end) + half_len;
for m = last_set_index + 1 : SignalLen
    M(:, m) = M(end:-1:1, 2*last_set_index - m + 1); % Repeat the pattern after the first beat
end

% test that all columns have at least one non-zero element:
% find(sum(M, 1) == 0)

data_prior_est = zeros(size(data));
data_posterior_est = zeros(size(data));
for ch = 1 : size(data, 1)
    x = data(ch, :); % the active channel

    stacked_beats = event_stacker(x, peak_indexes, event_width);
    [ECG_mean, ECG_var_mn, ECG_median, ECG_var_md] = robust_weighted_average(stacked_beats);
    sample_indexes = (0 : size(stacked_beats, 2) - 1) - (size(stacked_beats, 2)/2);

    switch params.BEAT_AVG_METHOD
        case 'MEAN'
            ECG_avg = ECG_mean;
            ECG_std = sqrt(ECG_var_mn);
        case 'MEDIAN'
            ECG_avg = ECG_median;
            ECG_std = sqrt(ECG_var_md);
    end




    % Under-development. For QRS and ventricular wave supression (useful for P-wave detection)
    % if isfield(params, 'avg_beat_shaping_window') && ~isempty(params.avg_beat_shaping_window)
    if isfield(params, 'Q_onset_lead') && ~isempty(params.Q_onset_lead) && ...
            isfield(params, 'T_offset_lag') && ~isempty(params.T_offset_lag)

        %{
        if(mod(params.avg_beat_shaping_window, 2) == 0)
            params.avg_beat_shaping_window = params.avg_beat_shaping_window + 1;
        end
        %}

        %{
        params.avg_beat_shaping_window = min(event_width, params.avg_beat_shaping_window);
        zero_padding_wings = (event_width - params.avg_beat_shaping_window) / 2;
        window = hamming(params.avg_beat_shaping_window);
        shape_window = [zeros(1, zero_padding_wings), window', zeros(1, zero_padding_wings)];
        %}

        %{
        nn_ = (1:length(ECG_avg));
        shape_window = exp(-(nn_ - ceil(length(ECG_avg)/2)).^2/(2*params.avg_beat_shaping_window^2));
        %}
        shape_window = SigmoidFunction(length(ECG_avg), ceil(length(ECG_avg)/2) - params.Q_onset_lead, 2.0) - ...
            SigmoidFunction(length(ECG_avg), ceil(length(ECG_avg)/2) + params.T_offset_lag, 0.5);

        % figure
        % plot(ECG_avg)
        % hold on
        % plot(shape_window*200)

        %         first_few_samples = 3;
        %         ECG_avg = 0        *      mean(ECG_avg(1:first_few_samples)) + ECG_avg .* shape_window;

        ECG_avg = ECG_avg .* shape_window;
    end
    % plot(ECG_avg)
    % grid


    switch params.NOISE_VAR_EST_METHOD
        case 'MIN' % Method 1: min std
            noise_std_est = min(ECG_std);
        case 'AVGLOWER' % Method 2: average of the smallest std
            if ~isfield(params, 'avg_bins')
                params.avg_bins = 3;
                warning(['Undefined numer of averaging bins; using ', num2str(params.avg_bins) ' as default']);
            end
            [ECG_std_up_sorted, ~] = sort(ECG_std);
            bns = 1 : params.avg_bins;
            noise_std_est = sqrt(sum(ECG_std_up_sorted(bns).^2 * (NumBeats - 1)) / (NumBeats * length(bns) - 1)); % Note: recovers the original variances and renormalizes by total (N-1) to obtain an unbiased estimator
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

    s_equiv_vars = zeros(1, SignalLen);
    for nn = 1 : SignalLen
        non_zeros_bins = find(M(:, nn));
        if isempty(non_zeros_bins)
            error('Empty column in mapping matrix');
        elseif length(non_zeros_bins) == 1 % only one point
            data_prior_est(ch, nn) = ECG_avg(non_zeros_bins);
            s_equiv_vars(nn) = ECG_intrinsic_var(non_zeros_bins);
        else % more than one point to average
            partial_prod_vars = zeros(1, length(non_zeros_bins));
            for bb = 1 : length(non_zeros_bins)
                non_zeros_excluding_current_ensemble = non_zeros_bins;
                non_zeros_excluding_current_ensemble(bb) = [];
                partial_prod_vars(bb) = prod(ECG_intrinsic_var(non_zeros_excluding_current_ensemble));
            end
            sum_partial_prod_vars = sum(partial_prod_vars);
            if sum_partial_prod_vars > 0
                data_prior_est(ch, nn) = sum(ECG_avg(non_zeros_bins(bb)) .* partial_prod_vars) / sum_partial_prod_vars;
                s_equiv_vars(nn) = prod(ECG_intrinsic_var(non_zeros_bins)) / sum_partial_prod_vars;
            else
                data_prior_est(ch, nn) = mean(ECG_avg(non_zeros_bins));
                s_equiv_vars(nn) = 0;
            end
        end

    end


    %     data_prior_est(ch, :) = (ECG_avg * M) ./ (ones(1, event_width) * M);




    %%%s_var = ECG_intrinsic_var * M_smoothed; % ECG beat variance estimate repeated over time

    % MAP estimate of each ECG sample assuming a Gaussian distribution for
    % the ECG samples and the noise (derived theoretically)
    data_posterior_est(ch, :) = (s_equiv_vars .* (x - params.n_mean) + n_var * data_prior_est(ch, :)) ./ (s_equiv_vars + n_var);

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
            errorbar(sample_indexes, ECG_avg, ECG_std/2);
            hold on
            plot(sample_indexes, ECG_avg);
            grid
            legend('Errorbar', 'Average ECG');
            title('Average ECG beat \pm STD');
            xlabel('Phase (rad)');
            ylabel('Amplitude');
            set(gca, 'fontsize', 16)

            figure
            plot(sample_indexes, ECG_std);
            hold on
            plot(sample_indexes, sqrt(ECG_intrinsic_var));
            plot(sample_indexes, sqrt(n_var)*ones(1, params.bins));
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
end

function y = SigmoidFunction(len, center, alpha)
n = (1 : len);
y = 1 ./ (1 + exp(-alpha*(n-center)));
end


