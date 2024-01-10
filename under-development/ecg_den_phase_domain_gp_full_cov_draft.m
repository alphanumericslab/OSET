function [data_posterior_est, data_prior_est, n_var] = ecg_den_phase_domain_gp_full_cov_draft(data, I_peaks, params)

if ~isfield(params, 'bins') % number of phase bins
    params.bins = round(max(diff(I_peaks))); % use median of RR-intervals in samples by default
end
if mod(params.bins, 2) == 0
    params.bins = params.bins + 1;
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

peaks = zeros(1, length(data));
peaks(I_peaks) = 1;

% Left padding the peaks sequence to add an additional peak on the far left
% (required for first beat correction)
first_rr_interval = I_peaks(2) - I_peaks(1);
r0_index = I_peaks(1) - first_rr_interval;
if r0_index < 1
    left_padding = 1 - r0_index;
    peaks_padded = cat(2, zeros(1, left_padding), peaks);
else
    left_padding = 0;
    peaks_padded = peaks;
end
% peaks_padded(I_peaks(1) + left_padding - first_rr_interval) = 1;
II_peaks = find(peaks_padded);

% Right padding the peaks sequence to add an additional peak on the far
% right (required for last beat correction)
last_rr_interval = II_peaks(end) - II_peaks(end-1);
r_last_index = II_peaks(end) + last_rr_interval;
if r_last_index > length(peaks_padded)
    right_padding = r_last_index - length(peaks_padded);
    peaks_padded = cat(2, peaks_padded, zeros(1, right_padding));
else
    right_padding = 0;
end
% II_peaks = find(peaks_padded);
% peaks_padded(I_peaks(end) + left_padding + last_rr_interval) = 1;
II_peaks = find(peaks_padded);


data_prior_est = zeros(size(data));
data_posterior_est = zeros(size(data));
data_padded = [zeros(size(data, 1), left_padding), data, zeros(size(data, 1), right_padding)];
for ch = 1 : size(data, 1)
    x_stacked = zeros(length(II_peaks), params.bins);
    peaks_mid_points = round((II_peaks(1:end-1) + II_peaks(2:end)) / 2);
    potential_first_mid_point = peaks_mid_points(1) - (peaks_mid_points(2) - peaks_mid_points(1));
    peaks_mid_points = cat(2, potential_first_mid_point, peaks_mid_points);

    potential_last_mid_point = peaks_mid_points(end) + (peaks_mid_points(end) - peaks_mid_points(end-1));
    peaks_mid_points = cat(2, peaks_mid_points, potential_last_mid_point);

    Phi = zeros(peaks_mid_points(1)-1, params.bins);
    phi_k = cell(1, length(II_peaks));
    for k = 1 : length(II_peaks)
        first_index = peaks_mid_points(k);
        last_index = peaks_mid_points(k+1) - 1;
        beat = data_padded(ch, first_index : last_index);
        theta_k = warping_transform_col([first_index, II_peaks(k), last_index] - first_index + 1, [1, ceil(params.bins/2), params.bins], params.interp_order);
        phi_k{k} = pinv(theta_k);
        % phi_k = (theta_k'*theta_k)\theta_k';
        Phi = cat(1, Phi, phi_k{k});
        x_stacked(k, :) = theta_k*beat';
    end
    Phi = cat(1, Phi, zeros(size(data_padded, 2) - size(Phi, 1), params.bins));









    Phi = Phi(left_padding + 1 : left_padding + size(data, 2), :);










    % if 0
    %     if left_padding > 0
    %         ECG_mean = mean(x_stacked(2 : length(II_peaks) - 1, :), 1);
    %         ECG_median = median(x_stacked(2 : length(II_peaks) - 1, :), 1);
    %         %     x_stacked_zero_mean = x_stacked - ECG_avg(ones(1, size(x_stacked, 1)), :);
    %         %         K_x_ph = x_stacked_zero_mean' * x_stacked_zero_mean / size(x_stacked_zero_mean, 1);
    %         %         K_x_ph = x_stacked(2:end-1, :)' * x_stacked(2:end-1, :) / (size(x_stacked, 1) - 2);
    %         K_x_ph = cov(x_stacked(2 : length(II_peaks) - 1, :)); % identical to the above lines if called as: cov(x_stacked, 1)
    %     else
    %         ECG_mean = mean(x_stacked(1 : length(II_peaks) - 1, :), 1);
    %         ECG_median = median(x_stacked(1 : length(II_peaks) - 1, :), 1);
    %         %     x_stacked_zero_mean = x_stacked - ECG_avg(ones(1, size(x_stacked, 1)), :);
    %         %         K_x_ph = x_stacked_zero_mean' * x_stacked_zero_mean / size(x_stacked_zero_mean, 1);
    %         %         K_x_ph = x_stacked(1:end-1, :)' * x_stacked(1:end-1, :) / (size(x_stacked, 1) - 1);
    %         K_x_ph = cov(x_stacked(1 : length(II_peaks) - 1, :)); % identical to the above lines if called as: cov(x_stacked, 1)
    %     end
    % end

    % peaks_temp = zeros(1, length(data_padded));
    % peaks_temp(II_peaks) = 1;
    % [phase, ~] = phase_calculator(peaks_temp); % phase calculation
    % [ECG_mean_0, ECG_std_0, meanphase_0, ECG_median_0, ECGSamplesPerBin_0] = avg_beat_calculator_phase_domain(data_padded, phase, params.bins, 0); % mean ECG extraction


    [ECG_mean, ~, ECG_median, ~] = robust_weighted_average(x_stacked);
    K_x_ph = cov(x_stacked, 1); % identical to the above lines if called as: cov(x_stacked, 1)

    switch params.BEAT_AVG_METHOD
        case 'MEAN'
            ECG_avg = ECG_mean;
        case 'MEDIAN'
            ECG_avg = ECG_median;
    end

    if isfield(params, 'nvar')
        n_var = params.nvar;
    else
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
        n_var = params.nvar_factor * noise_std_est^2; % noise variance estimate
    end

    n_var = 1.0 * n_var;




    % smooth the phase-time mapping matrix in time and phase
    switch params.SMOOTH_COV
        case 'BYPASS'

        case 'MA'
            if ~isfield(params, 'smooth_cov_wlen')
                params.smooth_cov_wlen = 3;
                warning(['Undefined kernel smoothing window length; using ', num2str(params.smooth_cov_wlen) ' as default']);
            end
            K_x_ph = filter2(ones(params.wlen_phase, params.wlen_time)/(params.wlen_phase * params.wlen_time), K_x_ph);
        case 'GAUSSIAN'
            if ~isfield(params, 'smooth_cov_gaussianstd')
                params.smooth_cov_gaussianstd = 3;
                warning(['Undefined gaussian smoothing window length; using ', num2str(params.smooth_cov_gaussianstd) ' as default']);
            end
            K_x_ph = imgaussfilt(K_x_ph, params.smooth_cov_gaussianstd);
    end
    K_x_ph = (K_x_ph + K_x_ph') / 2;


    % K_s_ph = K_x_ph - n_var * eye(size(K_x_ph));
    [V, d] = eig(K_x_ph, "vector");
    d = d - n_var;
    d(d < 0) = 0;
    K_s_ph = V * diag(d) * V';






    % K_s_ph = diag(diag(K_s_ph))







    METHOD = 'PHASE_DOMAIN_GP_IMPL'

    switch METHOD
        case '-1'
            K_s = Phi * K_s_ph * Phi';


            % K_s = imgaussfilt(K_s, params.smooth_cov_gaussianstd);


            K_s = (K_s + K_s') / 2;
            [U, d] = eig(K_s, "vector");
            % d(d < 0) = 0;
            rho = d / n_var;
            A = U * diag(rho./(1 + rho)) * U';
            I_A = U * diag(1./(1 + rho)) * U';
            % A = zeros(length(data_padded(ch, :)));
            % I_A = A;
            % for ii = 1 : length(rho)
            %     uut = U(:, ii)*U(:, ii)';
            %     A = A +  rho(ii)/(1 + rho(ii)) * uut;
            %     I_A = I_A +  1/(1 + rho(ii)) * uut;
            % end
            A = (A + A') / 2.0; % force symmetric if otherwise, due to computational errors
            I_A = (I_A + I_A') / 2.0; % force symmetric if otherwise, due to computational errors
            data_prior_est_padded = Phi * ECG_avg';
            % data_post_est_padded = I_A * data_prior_est_padded + A * data_padded(ch, :)';
            data_post_est_padded = I_A * data_prior_est_padded + A * data(ch, :)';
        case '0'
            K_s = Phi * K_s_ph * Phi';
            n_rank = length(I_peaks);
            [U, D] = eigs(K_s, n_rank);
            D = abs(diag(D));
            U = real(U);
            I_non_zero = find(10*log10(D) > -60.0);
            A = zeros(length(data_padded(ch, :)));
            for ii = 1 : length(I_non_zero)
                rho_i = D(I_non_zero(ii)) / n_var;
                if rho_i < 1e6 %isfinite(rho_i)
                    A = A +  rho_i/(1 + rho_i) * U(:, I_non_zero(ii))*U(:, I_non_zero(ii))';
                else
                    A = A +  1.0 * U(:, I_non_zero(ii))*U(:, I_non_zero(ii))';
                end
            end
            A = (A + A') / 2.0; % force symmetric if otherwise, due to computational errors
            data_prior_est_padded = Phi * ECG_avg';
            % data_post_est_padded = (eye(size(data_padded, 2)) - A) * data_prior_est_padded + A * data_padded(ch, :)';
            data_post_est_padded = (eye(size(data_padded, 2)) - A) * data_prior_est_padded + A * data(ch, :)';
        case '1'
            K_s = Phi * K_s_ph * Phi';
            K_x = K_s + n_var*eye(size(data_padded, 2));
            A = K_s * pinv(K_x);
            A = (A + A') / 2.0; % force symmetric if otherwise, due to computational errors
            data_prior_est_padded = Phi * ECG_avg';
            data_post_est_padded = (eye(size(data_padded, 2)) - A) * data_prior_est_padded + A * data_padded(ch, :)';
        case '2' % first method
            rho = K_s_ph / n_var; % SNR matrix
            Gamma = Phi * pinv(pinv(rho) + Phi'*Phi) * Phi';
            Gamma = (Gamma + Gamma') / 2.0; % force symmetric if otherwise, due to computational errors
            A = (Phi*rho*Phi')*(eye(size(data_padded, 2)) - Gamma);
            A = (A + A') / 2.0; % force symmetric if otherwise, due to computational errors

            data_prior_est_padded = Phi * ECG_avg';
            % data_post_est_padded = data_prior_est_padded + A*(data_padded(ch, :)' - data_prior_est_padded);
            data_post_est_padded = data_prior_est_padded + A*(data(ch, :)' - data_prior_est_padded);
        case '3' % second method
            rho = K_s_ph / n_var; % SNR matrix
            [V, d] = eig(rho, "vector");
            d(d < 0) = 0;
            S = V * diag(sqrt(d)) * V'; % symmetric semi-positive definite square root of rho
            S = (S + S') / 2;
            Upsilon = Phi * S;
            [U, sigma, ~] = svd(Upsilon);
            sigma = diag(sigma);
            % sigma(sigma < 0) = 0;
            A = zeros(length(data_padded(ch, :)));
            I_A = A;
            for ii = 1 : length(sigma)
                A = A + sigma(ii)^2 / (1 + sigma(ii)^2) * U(:, ii)*U(:, ii)';
                I_A = I_A + 1 / (1 + sigma(ii)^2) * U(:, ii)*U(:, ii)';
            end
            data_prior_est_padded = Phi * ECG_avg';
            % data_post_est_padded = (eye(size(data_padded, 2)) - A) * data_prior_est_padded + A * data_padded(ch, :)';
            % data_post_est_padded = I_A * data_prior_est_padded + A * data_padded(ch, :)';
            data_post_est_padded = I_A * data_prior_est_padded + A * data(ch, :)';
            norm_I_minus_A = norm(I_A)
            norm_A = norm(A)
        case 'PHASE_DOMAIN_GP_IMPL'
            data_prior_est_padded = [];%zeros(peaks_mid_points(1)-1, 1);
            data_post_est_padded = [];%zeros(peaks_mid_points(1)-1, 1);
            Ks_Kx_inv = K_s_ph*pinv(K_x_ph);
            Ks_Kx_inv = (Ks_Kx_inv + Ks_Kx_inv') / 2;
            for k = 1 : length(II_peaks)
                data_prior_est_padded = cat(1, data_prior_est_padded, phi_k{k}*ECG_avg');
                s_hat = Ks_Kx_inv * (x_stacked(k, :) - ECG_avg)' + ECG_avg';
                data_post_est_padded = cat(1, data_post_est_padded, phi_k{k}*s_hat);
            end
            data_prior_est_padded = cat(1, zeros(peaks_mid_points(1)-1, 1), data_prior_est_padded);
            data_post_est_padded = cat(1, zeros(peaks_mid_points(1)-1, 1), data_post_est_padded);
            
            size(data, 2)
            left_padding
            size(data_prior_est_padded)
            data_prior_est_padded = data_prior_est_padded(left_padding + 1 : left_padding + size(data, 2))';
            data_post_est_padded = data_post_est_padded(left_padding + 1 : left_padding + size(data, 2))';
    end

    data_prior_est(ch, :) = data_prior_est_padded;%(left_padding + 1 : left_padding + size(data, 2))';
    data_posterior_est(ch, :) = data_post_est_padded;%(left_padding + 1 : left_padding + size(data, 2))';

    lgnd = {};
    figure
    hold on
    % plot(data_padded(ch, :), 'linewidth', 2); lgnd = cat(2, lgnd, 'x');
    plot(data(ch, :), 'linewidth', 2); lgnd = cat(2, lgnd, 'x');
    plot(data_prior_est_padded, 'linewidth', 2); lgnd = cat(2, lgnd, '\mu_s');
    plot(data_post_est_padded, 'linewidth', 2); lgnd = cat(2, lgnd, '(I - A)*\mu_s + A*x');
    % plot((eye(size(data_padded, 2)) - A) * data_prior_est_padded, 'linewidth', 1); lgnd = cat(2, lgnd, '(I - A)*\mu');
    % plot((eye(size(data, 2)) - A) * data_prior_est_padded, 'linewidth', 1); lgnd = cat(2, lgnd, '(I - A)*\mu');
    % plot(A * data(ch, :)', 'linewidth', 1); lgnd = cat(2, lgnd, 'A*x');
    legend(lgnd);
    grid
    set(gca, 'fontsize', 16)

end

if 1
    figure
    plot(x_stacked')
    grid
    nn_padded = (0 : length(data_padded(ch, :)) - 1);
    figure
    plot(nn_padded, data_padded(ch, :));
    hold on
    plot(nn_padded(II_peaks), data_padded(ch, II_peaks), 'ro', 'markersize', 18);
    plot(nn_padded(peaks_mid_points), data_padded(ch, peaks_mid_points), 'gx', 'markersize', 18);
    grid

    lgnd = {};
    nn = (0 : length(data(ch, :)) - 1);
    figure
    plot(nn, data(ch, :)); lgnd = cat(2, lgnd, 'data');
    hold on
    plot(nn, data_prior_est(ch, :)); lgnd = cat(2, lgnd, 'data_prior_est');
    plot(nn, data_posterior_est(ch, :)); lgnd = cat(2, lgnd, 'data_posterior_est');
    grid
    legend(lgnd, 'interpreter', 'none');

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

% function y = LinearWarp(x, L)
% if isvector(x)
%     M = length(x);
%     tx = (0 : M - 1) / (M - 1);
%     ty = (0 : L - 1) / (L - 1);
%     y = interp1(tx, x, ty);
% elseif ismatrix(x)
%     M1 = size(x, 1);
%     M2 = size(x, 2);
%     [t1, t2] = meshgrid((0 : M2 - 1) / (M2 - 1), (0 : M1 - 1) / (M1 - 1));
%     [ty1, ty2] = meshgrid((0 : L(2) - 1) / (L(2) - 1), (0 : L(1) - 1) / (L(1) - 1));
%     y = interp2(t1, t2, x, ty1, ty2);
% else
%     error('First input should be either a vector or a matrix');
% end
% end
