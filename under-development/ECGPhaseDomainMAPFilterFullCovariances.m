function [data_posterior_est, data_prior_est, n_var] = ECGPhaseDomainMAPFilterFullCovariances(data, peaks, params)

% % % [phase, ~] = PhaseCalculation(peaks); % phase calculation

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


%
I_peaks = find(peaks);

if 0
    peaks_mid_points = round((I_peaks(1:end-1) + I_peaks(2:end)) / 2);
    potential_first_mid_point = peaks_mid_points(1) - (peaks_mid_points(2) - peaks_mid_points(1));
    if potential_first_mid_point >= 1
        peaks_mid_points = cat(2, potential_first_mid_point, peaks_mid_points);
    end
    potential_last_mid_point = peaks_mid_points(end) + (peaks_mid_points(end) - peaks_mid_points(end-1));
    if potential_last_mid_point <= length(peaks)
        peaks_mid_points = cat(2, peaks_mid_points, potential_last_mid_point);
    end

    I_peaks = peaks_mid_points;
    peaks = zeros(1, length(peaks));
    peaks(peaks_mid_points) = 1;
end

% if 0
    % Left padding the peaks sequence to add an additional peak on the far left
    % (required for first beat correction)
    if I_peaks(1) > 1
        rr0 = I_peaks(2) - I_peaks(1);
        r0_index = I_peaks(1) - rr0;
        if r0_index < 1
            left_padding = 1 - r0_index;
            peaks_padded = [1, zeros(1, left_padding - 1), peaks];
        else
            left_padding = 1;
            peaks_padded = [1, peaks];
        end
    else
        left_padding = 0;
        peaks_padded = peaks;
    end
    II_peaks = find(peaks_padded);

    % Right padding the peaks sequence to add an additional peak on the far
    % right (required for last beat correction)
    if II_peaks(end) < length(peaks_padded)
        rr_last = II_peaks(end) - II_peaks(end-1);
        r_last_index = II_peaks(end) + rr_last;
        if r_last_index > length(peaks_padded)
            right_padding = r_last_index - length(peaks_padded);
            peaks_padded = cat(2, peaks_padded, [zeros(1, right_padding - 1), 1]);
        else
            right_padding = 1;
            peaks_padded = cat(2, peaks_padded, 1);
        end
    else
        right_padding = 0;
        % peaks_padded = peaks_padded;
    end
    II_peaks = find(peaks_padded);
% end

% left_padding = 0;
% right_padding = 0;
% II_peaks = find(peaks);





data_prior_est = zeros(size(data));
data_posterior_est = zeros(size(data));
data_padded = [zeros(size(data, 1), left_padding), data, zeros(size(data, 1), right_padding)];
for ch = 1 : size(data, 1)
    x_stacked = zeros(length(II_peaks), params.bins);
    for k = 2 : length(II_peaks)
        previous_peak = II_peaks(k - 1);
        current_peak = II_peaks(k);
        beat = data_padded(ch, previous_peak : current_peak);
        x_phase_warped = LinearWarp(beat, params.bins + 1);
        x_stacked(k - 1, :) = x_phase_warped(1 : end - 1);
    end


    if left_padding > 0
        ECG_mean = mean(x_stacked(2 : length(II_peaks) - 1, :), 1);
        ECG_median = median(x_stacked(2 : length(II_peaks) - 1, :), 1);
        %     x_stacked_zero_mean = x_stacked - ECG_avg(ones(1, size(x_stacked, 1)), :);
        %         K_x_ph = x_stacked_zero_mean' * x_stacked_zero_mean / size(x_stacked_zero_mean, 1);
        %         K_x_ph = x_stacked(2:end-1, :)' * x_stacked(2:end-1, :) / (size(x_stacked, 1) - 2);
        K_x_ph = cov(x_stacked(2 : length(II_peaks) - 1, :)); % identical to the above lines if called as: cov(x_stacked, 1)
    else
        ECG_mean = mean(x_stacked(1 : length(II_peaks) - 1, :), 1);
        ECG_median = median(x_stacked(1 : length(II_peaks) - 1, :), 1);
        %     x_stacked_zero_mean = x_stacked - ECG_avg(ones(1, size(x_stacked, 1)), :);
        %         K_x_ph = x_stacked_zero_mean' * x_stacked_zero_mean / size(x_stacked_zero_mean, 1);
        %         K_x_ph = x_stacked(1:end-1, :)' * x_stacked(1:end-1, :) / (size(x_stacked, 1) - 1);
        K_x_ph = cov(x_stacked(1 : length(II_peaks) - 1, :)); % identical to the above lines if called as: cov(x_stacked, 1)
    end
    switch params.BEAT_AVG_METHOD
        case 'MEAN'
            ECG_avg = ECG_mean;
        case 'MEDIAN'
            ECG_avg = ECG_median;
    end

    ECG_std = sqrt(diag(K_x_ph))';

    switch params.NOISE_VAR_EST_METHOD
        case 'MIN' % Method 0: min std
            noise_std_est = min(ECG_std);
            %         case 'MINEIG' % Method 1: minimum eigenvalue
        case 'AVGLOWER' % Method 2: average of the smallest std
            if ~isfield(params, 'avg_bins')
                params.avg_bins = 3;
                warning(['Undefined numer of averaging bins; using ', num2str(params.avg_bins) ' as default']);
            end
            [ECG_std_up_sorted, std_sorted_indexes] = sort(ECG_std);
            ECGSamplesPerBin_sorted = ECG_std_up_sorted(std_sorted_indexes);
            bns = 1 : params.avg_bins;
            noise_std_est = sqrt(sum(ECG_std_up_sorted(bns).^2 .* (ECGSamplesPerBin_sorted(bns) - 1)) / (sum(ECGSamplesPerBin_sorted(bns)) - 1)); % Note: recovers the original variances and renormalizes by total (N-1) to obtain an unbiased estimator
            %             noise_std_est = min(ECG_std);
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




    n_var = 1.0 * n_var




    K_s_ph = K_x_ph - n_var * eye(size(K_x_ph, 1));

    %     K_s = zeros(length(data(ch, :)));

    K_s_padded = zeros(length(peaks_padded));
    data_prior_est_padded = zeros(1, length(peaks_padded));
    for k = 2 : length(II_peaks)
        current_peak = II_peaks(k);
        previous_peak = II_peaks(k - 1);
        beat = LinearWarp(ECG_avg, current_peak - previous_peak);
        data_prior_est_padded(previous_peak : current_peak - 1) = beat;
        for m = 2 : length(II_peaks)
            current_peak2 = II_peaks(m);
            previous_peak2 = II_peaks(m - 1);
            KK = LinearWarp(K_s_ph, [current_peak - previous_peak, current_peak2 - previous_peak2]);
            K_s_padded(previous_peak : current_peak - 1, previous_peak2 : current_peak2 - 1) = KK;
        end

    end
    data_prior_est(ch, :) = data_prior_est_padded(left_padding + 1 : end - right_padding);
    K_s = K_s_padded(left_padding + 1 : end - right_padding, left_padding + 1 : end - right_padding);

    %     mask = diag([zeros(1, 80), ones(1, size(K_s, 1) - 160), zeros(1, 80)]);
    %     mask = imgaussfilt(mask, 250.0);
    %     mask = mask./max(abs(mask(:)));
    %
    %     % mask = max(mask(:)) - mask;
    %     K_s = K_s .* mask;

    if 0
        [U,S,V] = svd(K_s);
        S(S < 0) = eps;
        K_s = U * S * V';
    end

    %     ww = 3;
    %     K_s = filter2(ones(ww, ww)/(ww * ww), K_s);
    %     K_s = imgaussfilt(K_s, 1.0);

    K_s = (K_s + K_s') / 2;

    K_x = K_s + n_var * eye(size(K_s, 1));
    P = (K_s / K_x);
    data_posterior_est(ch, :) = data_prior_est(ch, :) + (data(ch, :) - data_prior_est(ch, :)) * P';
end

if 1
    figure
    plot(data)
    hold on
    plot(I_peaks, data(I_peaks), 'ro', 'markersize', 18)
    %     plot(data_posterior_est(ch, :))
    plot(data_prior_est(ch, :))
    %     plot(9*sqrt(abs((diag(K_s)))), 'k')
    %     legend('raw', 'peaks', 'posterior', 'prior', 'K_s');
    legend('raw', 'peaks', 'prior');
    % legend('raw', 'peaks', 'posterior', 'prior');
    grid

    figure
    plot(x_stacked')
    hold on
    plot(sqrt(diag(K_x_ph)), 'k', 'linewidth', 3)
    grid

    % figure
    % plot(P);
    % grid
end
% for ch = 1 : size(data, 1)
%     x = data(ch, :); % the active channel
%
%     % construct a stacked matrix representation of the ECG beats in the phased omain
%     x_stacked = zeros(length(end_of_beat_indexes), params.bins);
%     for ll = 1 : length(end_of_beat_indexes)
%         if ll > 1
%             beat_start = end_of_beat_indexes(ll-1) + 1;
%         else
%             beat_start = 1;
%         end
%         beat_end = end_of_beat_indexes(ll);
%         x_stacked(ll, :) = LinearWarp(x(beat_start : beat_end), params.bins);%x(beat_start : beat_end) * (diag(max(1, sum(M(:, beat_start : beat_end), 2))) \ M(:, beat_start : beat_end))';
%         %         x_stacked(ll, :) = x(beat_start : beat_end) * M(:, beat_start : beat_end)';
%     end
%
%     mu_s = [];
%     for ll = 1 : length(end_of_beat_indexes)
%         if ll > 1
%             beat_start = end_of_beat_indexes(ll-1) + 1;
%         else
%             beat_start = 1;
%         end
%         beat_end = end_of_beat_indexes(ll);
%         %             mu_s = cat(2, mu_s, LinearWarp(ECG_mean, beat_end - beat_start + 1));
%         data_prior_est(ch, beat_start : beat_end) = LinearWarp(ECG_mean, beat_end - beat_start + 1);
%     end
%
%     % disp(['nvar estimate = ' num2str(sqrt(n_var))])
%     %         y = LinearWarp(ECG_mean, 270);
%     %         z = LinearWarp(y, length(ECG_mean));
%     %
%     [~, ~, meanphase, ~, ECGSamplesPerBin] = MeanECGExtraction(x, phase, params.bins, 1); % mean ECG extraction
%
%     ECG_intrinsic_var = max(0, ECG_std.^2 - n_var); % average beat variance estimate
%
%     % smooth the phase-time mapping matrix in time and phase
%     switch params.SMOOTH_PHASE
%         case 'BYPASS'
%             M_smoothed = M;
%         case 'MA'
%             if ~isfield(params, 'wlen_time')
%                 params.wlen_time = 3;
%                 warning(['Undefined time smoothing window length; using ', num2str(params.wlen_time) ' as default']);
%             end
%             if ~isfield(params, 'wlen_phase')
%                 params.wlen_phase = 3;
%                 warning(['Undefined phase smoothing window length; using ', num2str(params.wlen_phase) ' as default']);
%             end
%             M_smoothed = filter2(ones(params.wlen_phase, params.wlen_time)/(params.wlen_phase * params.wlen_time), M);
%         case 'GAUSSIAN'
%             if ~isfield(params, 'gaussianstd')
%                 params.gaussianstd = 3;
%                 warning(['Undefined gaussian smoothing window length; using ', num2str(params.gaussianstd) ' as default']);
%             end
%             M_smoothed = imgaussfilt(M, params.gaussianstd);
%     end
%
%     % reconstruct the ECG
%     if 0 % under dev.
%         K_s_ph = K_x_ph - n_var * eye(params.bins);
%         rho = K_s_ph / n_var;
%         sig_len = length(x);
%         T = M_smoothed / length(end_of_beat_indexes);
%         %     K_s_tm = M' * K_s_ph * M;
%         %     K_x_tm = M' * K_s_ph * M + n_var * eye(size(M, 2));
%         data_prior_est(ch, :) = (0.5*length(end_of_beat_indexes)) * (T'*T) * x';%ECG_avg * T; % prior estimate based on average beat repetition over time
%         % MAP estimate of each ECG sample assuming a Gaussian distribution for
%         % the ECG samples and the noise (derived theoretically)
%         %     data_posterior_est(ch, :) = ((T' * rho * T) * (eye(sig_len) - T' * pinv(pinv(rho) + (T * T')) * T) * (eye(sig_len) - (T'*T)) + (T'*T)) * x';
%         data_posterior_est(ch, :) = (T' * rho * T) * (eye(sig_len) - T' * pinv(pinv(rho) + (T * T')) * T) * (x' - data_prior_est(ch, :)') + data_prior_est(ch, :)';
%     else % under dev.
%         K_s_ph = K_x_ph - n_var * eye(params.bins);
%
%         % make matrix SPD
%         [U,S,V] = svd(K_s_ph);
%         S(S < 0) = 0;
%         K_s_ph = U * S * V';
%
%         rho = K_s_ph / n_var;
%         sig_len = length(x);
%         T = M_smoothed;
%         %     K_s_tm = M' * K_s_ph * M;
%         %     K_x_tm = M' * K_s_ph * M + n_var * eye(s ize(M, 2));
%         %T2 = diag(sum(T, 2)) \ T;
%         T2 = T / diag(sum(T, 1));
%
%         %         data_prior_est(ch, :) = T' * (T2 * x');%ECG_avg * T; % prior estimate based on average beat repetition over time
%         data_prior_est(ch, :) = ECG_avg * T2;%ECG_avg * T; % prior estimate based on average beat repetition over time
%
%         % MAP estimate of each ECG sample assuming a Gaussian distribution for
%         % the ECG samples and the noise (derived theoretically)
%         %     data_posterior_est(ch, :) = ((T' * rho * T) * (eye(sig_len) - T' * pinv(pinv(rho) + (T * T')) * T) * (eye(sig_len) - (T'*T)) + (T'*T)) * x';
%
%         %         data_posterior_est(ch, :) = (T2' * rho * T2) * (eye(sig_len) - T2' * pinv(pinv(rho) + (T2 * T2')) * T2) * (x' - data_prior_est(ch, :)') + data_prior_est(ch, :)';
%
%         Delta = x' - data_prior_est(ch, :)';
%         term1 = T2 * Delta;
%         term2 = T2' * rho;
%         term3 = T2 * T2';
%         term3_ = tril(triu(term3, -round(params.bins/2)), round(params.bins/2)); % removes the phase glithces at the corners of term3
%         data_posterior_est(ch, :) = (term2 * term1) - term2 * term3_ * pinv(pinv(rho) + term3_) * term1 + data_prior_est(ch, :)';
%
%         if 0
%             figure
%             plot(x)
%             hold on
%             plot(data_prior_est(ch, :));
%             plot(data_posterior_est(ch, :));
%             grid
%         end
%     end
%
%     %     plot(x)
%     %     hold on
%     %     plot(data_prior_est(ch, :))
%     %     plot(data_posterior_est(ch, :))
%     %     grid
%
%     if isfield(params, 'plotresults')
%         if params.plotresults == true
%             % if ch == 1
%             %     figure
%             %     plot(phase)
%             %     hold on
%             %     plot(pphase)
%             %     grid
%             %     legend('phase', 'pphase');
%             % end
%
%             figure
%             errorbar(meanphase, ECG_avg, ECG_std/2);
%             hold on
%             plot(meanphase, ECG_avg);
%             grid
%             legend('Errorbar', 'Average ECG');
%             title('Average ECG beat \pm STD');
%             xlabel('Phase (rad)');
%             ylabel('Amplitude');
%             set(gca, 'fontsize', 16)
%
%             figure
%             plot(meanphase, ECG_std);
%             hold on
%             plot(meanphase, sqrt(ECG_intrinsic_var));
%             plot(meanphase, sqrt(n_var)*ones(1, params.bins));
%             grid
%             legend('Signal+Noise STD', 'Signal STD', 'Noise STD');
%             title('Signal and noise STDs');
%             xlabel('Phase (rad)');
%             ylabel('Amplitude');
%             set(gca, 'fontsize', 16)
%
%             figure
%             plot(sort(ECG_std));
%             grid
%             title('Sample STDs across the average beat (sorted)');
%             xlabel('Sorted bin index');
%             ylabel('STDd');
%             set(gca, 'fontsize', 16)
%
%             figure
%             plot(x)
%             hold on
%             plot(x, 'linewidth', 2)
%             plot(data_prior_est(ch, :), 'linewidth', 2)
%             plot(data_posterior_est(ch, :), 'linewidth', 2)
%             grid
%             legend('Noisy ECG', 'Prior ECG estimate', 'Posterior ECG estimate')
%             xlabel('time (samples)');
%             ylabel('Amplitude');
%             title('Filtering results');
%             set(gca, 'fontsize', 16)
%         end
%     end
% end
end

function y = LinearWarp(x, L)
if isvector(x)
    M = length(x);
    tx = (0 : M - 1) / (M - 1);
    ty = (0 : L - 1) / (L - 1);
    y = interp1(tx, x, ty);
elseif ismatrix(x)
    M1 = size(x, 1);
    M2 = size(x, 2);
    [t1, t2] = meshgrid((0 : M2 - 1) / (M2 - 1), (0 : M1 - 1) / (M1 - 1));
    [ty1, ty2] = meshgrid((0 : L(2) - 1) / (L(2) - 1), (0 : L(1) - 1) / (L(1) - 1));
    y = interp2(t1, t2, x, ty1, ty2);
else
    error('First input should be either a vector or a matrix');
end
end
