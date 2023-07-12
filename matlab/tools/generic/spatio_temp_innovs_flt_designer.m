function [h, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = spatio_temp_innovs_flt_designer(template_records, varargin)
% spatio_temp_innovs_flt_designer - Designs an innovations filter from ensembles of multichannel random processes.
%
% Syntax: [h, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = spatio_temp_innovs_flt_designer(template_records, params)
%
% Inputs:
%   template_records: A cell array or a single matrix used for spectral estimation. Each element of the cell array is considered as an ensemble
%                     of the multichannel random process of interest.
%   params:           A structure with the following fields
%         params.spatial_filter_type:            'BY_PASS', 'PCA', or 'ICA' - Spatial filter applied before innovations filter design
%         params.normalize_records:              true/false - Normalize channels or not
%         params.fs:                             Sampling frequency in Hz
%         params.keep_mean:                      true/false - Keep or remove the channel-wise means
%         params.spectral_len:                   Number of spectral estimation bins (512 by default)
%         params.filter_len:                     Innovations filter length (512 by default). Best practice to set filter_len = spectral_len
%                                                (but may result in long impulse responses)
%         params.smooth_spectrum:                true/false - Smooth the estimated spectrum before filter design or not
%         params.lambda:                         Tikhonov regularization parameter used for spectral smoothing (10000.0 by default)
%         params.spectral_averaging_method:       Spectral averaging method 'MEAN', 'MEDIAN', 'MAX', 'MIN', or 'MAX_MIN_AVG' (Default: 'MEDIAN')
%         params.innovation_filter_type:          'LINEAR_PHASE' or 'MIN_PHASE' (Default: 'LINEAR_PHASE'). In 'MIN_PHASE' mode, params.filter_len
%                                                should be odd. Note: 'MIN_PHASE' is slow for long filter lengths. Not recommended unless a
%                                                'minimum-phase' filter is specifically required.
%         params.plot_results:                   true/false - Plot the results or not
%
% Outputs:
%   h:                Cell array containing the innovations filters
%   A:                Transformation matrix used for spatial filtering
%   S_mean:           Mean spectral estimation
%   S_median:         Median spectral estimation
%   S_max:            Maximum spectral estimation
%   S_min:            Minimum spectral estimation
%   S_max_min_avg:    Average of maximum and minimum spectral estimations
%
% Algorithm: The spectrum of the multichannel data is estimated per channel using the Welch spectral estimation algorithm with a Hamming window.
%            An innovations filter is designed by factorizing the estimated spectra into linear or minimum-phase factors.
%
%   Revision History:
%       2023: First release
% 
%   Reza Sameni, 2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Load the default parameters
default_params = DefaultParameters();
if nargin > 1
    params = varargin{1};
else
    params = default_params;
end

incomplete_parameter_list = false;
% Check and assign parameter values
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

% Convert into cell format if input is a single matrix
if isvector(template_records)
    template_records = template_records(:)'; % Channels are row-wise
end
if ~iscell(template_records)
    template_records = {template_records};
end

N_channels = size(template_records{1}, 1); % Number of channels
for k = 2:length(template_records)
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

% Check and correct filter length for minimum-phase mode
if mod(params.filter_len, 2) == 0 && isequal(params.innovation_filter_type, 'MIN_PHASE')
    params.filter_len = params.filter_len + 1;
    warning(['Filter length needs to be odd in minimum-phase mode. New filter size: ', num2str(params.filter_len)]);
end

P = zeros(length(template_records), N_channels, params.spectral_len);
for k = 1:length(template_records)
    x = template_records{k};
    % Perform PCA, ICA, or bypass spatial filter
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

    % Bypass or apply record-wise/segment-wise normalization
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

% Remove channel-wise mean if keep_mean is false
if ~params.keep_mean
    S_mean(:, 1) = S_global_min;
    S_median(:, 1) = S_global_min;
    S_max(:, 1) = S_global_min;
    S_min(:, 1) = S_global_min;
    S_max_min_avg(:, 1) = S_global_min;
end

% Select the spectral averaging method
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

% Smooth the spectrum if required
if params.smooth_spectrum
    S_smoothed = zeros(size(S));
    for m = 1:N_channels
        S_smoothed(m, :) = tikhonov_regularization(S(m, :), 2, params.lambda);
        I = S_smoothed(m, :) < 0;
        S_smoothed(m, I) = 0;
        S_smoothed(m, :) = (S_smoothed(m, :) + S_smoothed(m, end:-1:1)) / 2;
    end
else
    S_smoothed = S;
end

h0 = real(ifft(sqrt(abs(S_smoothed)), params.filter_len, 2));
h0_zero = fftshift(h0, 2);
h0_zero = (h0_zero + h0_zero(:, end:-1:1)) / 2;

h = cell(1, N_channels);
switch params.innovation_filter_type
    case 'LINEAR_PHASE'
        for ch = 1:N_channels
            h(ch) = {h0_zero(ch, :)};
        end
    case 'MIN_PHASE'
        for ch = 1:N_channels
            r = conv(h0_zero(ch, :), h0_zero(ch, end:-1:1));
            h(ch) = {firminphase(r)};
        end
    otherwise
        error('Undefined innovations filter type');
end

% Plot the results if required
if params.plot_results
    ff = params.fs * (0:params.spectral_len-1) / params.spectral_len;
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
    xlabel('Frequency (Hz)');
    ylabel('Magnitude response (dB)');
    legend(lgnd);
    grid
    title('Estimated spectra using different methods');
    set(gca, 'fontsize', 18)

    for ch = 1:N_channels
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
        xlabel('Frequency (Hz)');
        ylabel('Magnitude response (dB)');
        sgtitle(['Innovations filter''s impulse response in Channel ' num2str(ch)], 'fontsize', 18);
        set(gca, 'fontsize', 18)
    end
end

end

% Function to generate default parameters
function params = DefaultParameters()
params.spatial_filter_type = 'BY_PASS'; % 'BY_PASS', 'PCA', or 'ICA'
params.normalize_records = true;
params.fs = 1.0;
params.keep_mean = true;
params.spectral_len = 512;
params.filter_len = 512;
params.lambda = 10000.0;
params.spectral_averaging_method = 'MEDIAN';
params.smooth_spectrum = true;
params.innovation_filter_type = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE'
params.plot_results = true;
end
