function [h, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = SpatioTemporalInnovationsFilterDesigner(template_records, varargin)

% [h, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = SpatioTemporalInnovationsFilterDesigner(template_records, params)
% Estimates an innovations filter from ensembles of multichannel random
% processes.
%
% Usage: innovations filters can be used to generate random processes that
% resemble a training dataset, when fed by input white noise.
%
% Inputs:
%   template_records: a cell array or a single matrix used for spectral
%       estimation. Each element of the cell array is considered as an ensemble
%       of the multichannel random process of interest.
%   params: a structure with the followoing fields
%         params.spatial_filter_type = 'BY_PASS', 'PCA' or 'ICA': spatial filter applied before innovations filter design
%         params.normalize_records = true/false: normalize channels or not
%         params.fs: sampling frequency in Hz
%         params.keep_mean = true/false: keep or remove the channel-wise means
%         params.spectral_len: number of spectral estimation bins (512 by default)
%         params.filter_len: innovations filter length (512 by default). Best practice to set filter_len = spectral_len (but may result in long impulse responses)
%         params.smooth_spectrum = true/false: smooth the estimated spectrum before filter design or not
%         params.lambda: Tikhonov regularization parameter used for spectral smoothing 10000.0
%         params.spectral_averaging_method: spectral averaging method 'MEAN', 'MEDIAN', 'MAX', 'MIN' or 'MAX_MIN_AVG'. Default: 'MEDIAN'
%         params.innovation_filter_type: 'LINEAR_PHASE' or 'MIN_PHASE'.
%         Default: 'LINEAR_PHASE'. In 'MIN_PHASE' mode params.filter_len should be odd. Note: 'MIN_PHASE' is slow for long filter lengths.
%           Not recommended unless a 'minimum-phase' filter is specifically required
%         params.plot_results = true/false: plot the results or not
%
% Algorithm: the spectrum of the multichannel data is estimated per
% channel, uisng the Welch spectral estimation algorithm with a Hamming
% window. An  innovations filter is designed by factorizing the estimated
% spectra into linear or minimum phase factors.
%
% The Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET
% Reza Sameni, Feb 2023

% Load the default parameters
default_params = DefaultParameters();
if nargin > 1
    params = varargin{1};
else
    params = default_params;
end

incomplete_parameter_list = false;
if ~isfield(params, 'spatial_filter_type')
    params.spatial_filter_type = default_params.spatial_filter_type;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'normalize_records')
    params.normalize_records = default_params.normalize_records;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'keep_mean')
    params.keep_mean = default_params.keep_mean;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'spectral_len')
    params.spectral_len = default_params.spectral_len;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'filter_len')
    params.filter_len = default_params.filter_len;
end
if ~isfield(params, 'smooth_spectrum')
    params.smooth_spectrum = default_params.smooth_spectrum;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'lambda') && params.smooth_spectrum == true
    params.lambda = default_params.lambda;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'spectral_averaging_method')
    params.spectral_averaging_method = default_params.spectral_averaging_method;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'innovation_filter_type')
    params.innovation_filter_type = default_params.innovation_filter_type;
    incomplete_parameter_list = true;
end
if ~isfield(params, 'plot_results')
    params.plot_results = default_params.plot_results;
end
if ~isfield(params, 'fs')
    params.fs = default_params.fs;
    if params.plot_results
        warning(['Sampling frequency not provided. Assuming fs = ', num2str(default_params.fs), '(Hz) by default in the plots.']);
    end
    incomplete_parameter_list = true;
end
if incomplete_parameter_list
    disp('Parameter list is incomplete. Using default parameters:');
    disp(params);
end



% Convert into cell format, if input is a single matrix
if isvector(template_records)
    template_records = template_records(:)'; % channels are row-wise
end

if ~iscell(template_records)
    template_records = {template_records};
end

N_channels = size(template_records{1}, 1); % number of channels
for k = 2 : length(template_records)
    if N_channels ~= size(template_records{k}, 1)
        error('All records should have the same number of channels.');
    end
end

% Preserve or remove the mean
if params.keep_mean
    mean_factor = 0;
else
    mean_factor = 1;
end

if mod(params.filter_len, 2) == 0 && isequal(params.innovation_filter_type, 'MIN_PHASE')
    params.filter_len = params.filter_len + 1;
    warning(['Filter length needs to be odd length in minimum-phase mode. New filter size: ', num2str(params.filter_len)]);
end

P = zeros(length(template_records), N_channels, params.spectral_len);
for k = 1 : length(template_records) % Sweep over all records

    x = template_records{k};
    % peform PCA, ICA, or bypass spatial filter
    switch params.spatial_filter_type
        case 'BY_PASS'
            x_decomposed = x;
            A = eye(N_channels);
        case 'PCA'
            [U, ~] = eig(cov(x'));
            x_decomposed = U' * x;
            A = U;
        case 'ICA'
            B = jadeR(x);
            x_decomposed = B * x;
            A = pinv(B);
    end

    % by-pass or apply record-wise/segment-wise normalization
    if params.normalize_records
        x_normalized = (x_decomposed - mean_factor * repmat(mean(x_decomposed, 2), 1, size(x_decomposed, 2))) ./ repmat(std(x_decomposed, [], 2), 1, size(x_decomposed, 2));
    else
        x_normalized = x_decomposed;
    end

    P(k, :, :) = pwelch(x_normalized', [], [], params.spectral_len, 'twosided')';
end

S_mean = shiftdim(mean(P, 1), 1);
S_median = shiftdim(median(P, 1), 1);
S_max = shiftdim(max(P, [], 1), 1);
S_min = shiftdim(min(P, [], 1), 1);
S_max_min_avg = (S_max + S_min) / 2;

% This is needed since spectral windowing can add DC
S_global_min = min(S_min(:)); % take this as the noise floor
if ~params.keep_mean
    S_mean(:, 1) = S_global_min;
    S_median(:, 1) = S_global_min;
    S_max(:, 1) = S_global_min;
    S_min(:, 1) = S_global_min;
    S_max_min_avg(:, 1) = S_global_min;
end


switch params.spectral_averaging_method
    case 'MEAN'
        S = S_mean;
    case 'MEDIAN'
        S = S_median;
    case 'MAX'
        S = S_max;
    case 'MIN'
        S = S_min;
    case 'MAX_MIN_AVG'
        S = S_max_min_avg;
end

if params.smooth_spectrum
    S_smoothed = zeros(size(S));
    for m = 1 : N_channels
        S_smoothed(m, :) = TikhonovRegularization(S(m, :), 2, params.lambda);
        I = S_smoothed(m, :) < 0;
        S_smoothed(m, I) = 0;
        S_smoothed(m, :) = (S_smoothed(m, :) + S_smoothed(m, end:-1:1))/2;
    end
else
    S_smoothed = S;
end

% S_complex = sqrt(S_smoothed) .* (ones(N_channels, 1) * exp(1j*2*pi*(0:params.spectral_len-1)));
h0 = real(ifft(sqrt(abs(S_smoothed)), params.filter_len, 2));
h0_zero = fftshift(h0, 2);
h0_zero = (h0_zero + h0_zero(:, end:-1:1)) / 2;
h = cell(1, N_channels);
switch params.innovation_filter_type
    case 'LINEAR_PHASE'
        for ch = 1 : N_channels
            h(ch) = {h0_zero(ch, :)};
        end
    case 'MIN_PHASE'
        for ch = 1 : N_channels
            r = conv(h0_zero(ch, :), h0_zero(ch, end:-1:1));
            h(ch) = {firminphase(r)};
        end
    otherwise
        error('Undefined innovations filter type');
end

% PLOT RESULTS
if params.plot_results
    ff = params.fs * (0:params.spectral_len-1)/params.spectral_len;
    scale = 1;
    ref = 0;

    lgnd = {};
    figure
    hold on
    plot(ff - params.fs/2, fftshift(10*log10(scale*S_max) - ref), 'linewidth', 3); lgnd = cat(1, lgnd, 'Max');
    plot(ff - params.fs/2, fftshift(10*log10(scale*S_min) - ref), 'linewidth', 3); lgnd = cat(1, lgnd, 'Min');
    plot(ff - params.fs/2, fftshift(10*log10(scale*S_mean) - ref), 'linewidth', 3); lgnd = cat(1, lgnd, 'Mean');
    plot(ff - params.fs/2, fftshift(10*log10(scale*S_median) - ref), 'linewidth', 3); lgnd = cat(1, lgnd, 'Medium');
    plot(ff - params.fs/2, fftshift(10*log10(scale*S_max_min_avg) - ref), 'linewidth', 3); lgnd = cat(1, lgnd, 'Min-Max Avg.');
    plot(ff - params.fs/2, fftshift(10*log10(scale*S_smoothed) - ref), 'linewidth', 3); lgnd = cat(1, lgnd, ['Smoothed ' params.spectral_averaging_method]);
    xlabel('frequency(Hz)');
    ylabel('Magnitude response (dB)');
    legend(lgnd);
    grid
    title('Estimated spectra using different methods');
    set(gca, 'fontsize', 18)

    for ch = 1 : N_channels
        [H, F] = freqz(h{ch}, 1, 'whole');

        figure
        subplot(211);
        stem(h{ch});
        grid
        xlabel('n');
        ylabel('Impulse response');
        set(gca, 'fontsize', 18)

        subplot(212)
        plot(ff, 10*log10(S_smoothed(ch, :)));
        hold on
        plot(params.fs * F/2/pi, 20*log10(abs(H)));
        grid
        legend('Data spectrum', 'Innovations filter''s squared magnitude response');
        xlabel('frequency(Hz)');
        ylabel('Magnitude response (dB)');
        sgtitle(['Innovations filter''s impulse response in Channel ' num2str(ch)], 'fontsize', 18);
        set(gca, 'fontsize', 18)
    end
end

end


function params = DefaultParameters()
params.spatial_filter_type = 'BY_PASS'; % 'BY_PASS', 'PCA' or 'ICA'
params.normalize_records = true;
params.fs = 1.0;
params.keep_mean = true;
params.spectral_len = 512;
params.filter_len = 512;
params.lambda = 10000.0;
params.spectral_averaging_method = 'MEDIAN';
params.smooth_spectrum = true;
params.innovation_filter_type = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE';
params.plot_results = true;
end
