function [erps_rwa, erps_mean, erps_rwmd, erps_var, stacked_targets] = erp_extractor_rwa(data, fs, target_indexes, erp_wlen, varargin)
% Event related potential extraction using Robust Weighted Averaging
%
% Reza Sameni (reza.sameni@gmail.com)
% Copyright June, 2019
%

% check if custom params have been provided
if nargin > 4 && ~isempty(varargin{1})
    params = varargin{1};
else
    params = struct;
end

if ~isfield(params, 'verbose') || isempty(params.verbose)
    params.verbose = true;
end
if params.verbose, disp('displaying all default parameters in verbose mode; set params.verbose=false for silent mode'), end

if ~isfield(params, 'plot_results') || isempty(params.plot_results)
    params.plot_results = false;
    if params.verbose, disp(['params.plot_results = ', params.plot_results]), end
end

if ~isfield(params, 'plots_per_figure') || isempty(params.plots_per_figure)
    params.plots_per_figure = 4;
    if params.verbose, disp(['params.plots_per_figure = ', num2str(params.plots_per_figure)]), end
end

if ~isfield(params, 'baseline_filter') || isempty(params.baseline_filter)
    params.baseline_filter = 'DC'; % 'NONE', 'DC', 'MDMN', 'SINGLE-ORDER-IIR'
    if params.verbose, disp(['params.baseline_filter = ', params.baseline_filter]), end
end

switch params.baseline_filter
    case 'MDMN'
        if ~isfield(params, 'baseline_wlen1') || isempty(params.baseline_wlen1)
            params.baseline_wlen1 = 0.3;
            if params.verbose, disp(['params.baseline_wlen1 = ', num2str(params.baseline_wlen1)]), end
        end

        if ~isfield(params, 'baseline_wlen2') || isempty(params.baseline_wlen2)
            params.baseline_wlen2 = 0.45;
            if params.verbose, disp(['params.baseline_wlen2 = ', num2str(params.baseline_wlen2)]), end
        end
        baseline_wlen1 = round(params.baseline_wlen1*fs); % convert into samples
        baseline_wlen2 = round(params.baseline_wlen2*fs); % convert into samples
    case 'SINGLE-ORDER-IIR'
        if ~isfield(params, 'hp_cuttoff') || isempty(params.hp_cuttoff)
            params.hp_cuttoff = 0.5;
            if params.verbose, disp(['params.hp_cuttoff = ', num2str(params.hp_cuttoff)]), end
        end
end

if ~isfield(params, 'enhancer_filter') || isempty(params.enhancer_filter)
    params.enhancer_filter = false;
    if params.verbose, disp(['params.enhancer_filter = ', num2str(params.enhancer_filter)]), end
end

if params.enhancer_filter
    if ~isfield(params, 'enhancer_filt_b') || isempty(params.enhancer_filt_b)
        error('params.enhancer_filt_b is required when params.enhancer_filter is true')
    end
    if ~isfield(params, 'enhancer_filt_a') || isempty(params.enhancer_filt_a)
        error('params.enhancer_filt_a is required when params.enhancer_filter is true')
    end
end

if ~isfield(params, 'enhancer_nsca') || isempty(params.enhancer_nsca)
    params.enhancer_nsca = false;
    if params.verbose, disp(['params.enhancer_nsca = ', params.enhancer_nsca]), end
end

if params.enhancer_nsca
    if ~isfield(params, 'nsca_max_components') || isempty(params.nsca_max_components)
        params.nsca_max_components = 1;
        if params.verbose, disp(['params.nsca_max_components = ', num2str(params.nsca_max_components)]), end
    end
end

if ~isfield(params, 'erp_vertical_offset_alignment') || isempty(params.erp_vertical_offset_alignment)
    params.erp_vertical_offset_alignment = 'NONE'; % 'NONE', 'ZERO-MEAN', 'ZERO-MEAN-FRACTION'
    if params.verbose, disp(['params.erp_vertical_offset_alignment = ', params.erp_vertical_offset_alignment]), end
end

if isequal(params.erp_vertical_offset_alignment , 'ZERO-MEAN-FRACTION')
    params.erp_vertical_offset_fraction = 0.1;
end

num_ch = size(data, 1); % number of channels
sig_len = size(data, 2); % time length

switch params.baseline_filter
    case 'NONE'
        baseline = zeros(size(data));
    case 'DC' % remove mean
        baseline = repmat(mean(data, 2), 1, size(data, 2));
    case 'MDMN' % remove baseline using two-step moving median and moving average filters
        bl = baseline_sliding_window(data, baseline_wlen1, 'md');
        baseline = baseline_sliding_window(bl, baseline_wlen2, 'mn');
    case 'SINGLE-ORDER-IIR' % single pole forward-backward IIR filter
        baseline = lp_filter_zero_phase(data, params.hp_cuttoff/fs);
end
data_detrended = data - baseline;

% denoise the signal if demanded
if params.enhancer_filter
    data_den = filtfilt(params.enhancer_filt_b, params.enhancer_filt_a, data_detrended')';
else
    data_den = data_detrended;
end

% apply nonstationary component analysis if demanded
if params.enhancer_nsca
    spikes = zeros(1, sig_len);
    spikes(target_indexes) = 1;
    pulse = filter(ones(1, erp_wlen), 1, spikes);
    pulse = find(pulse > 0);
    [y, ~, A] = nonstationary_component_analysis(data_den, pulse, 1:sig_len);
    data_den = A(:, 1:params.nsca_max_components) * y(1:params.nsca_max_components, :);
end

% stack all ERP across all channels
stacked_targets = zeros(num_ch, length(target_indexes), erp_wlen); % matrix of stacked ERPs
for tgt_indx = 1 : length(target_indexes)
    start = target_indexes(tgt_indx);
    segment = data_den(:, start : min(start + erp_wlen - 1, sig_len));
    switch params.erp_vertical_offset_alignment
        case 'NONE'
            mn = 0;
        case 'ZERO-MEAN'
            mn = mean(segment, 2);
        case 'ZERO-MEAN-FRACTION'
            mn = mean(segment(:, 1:ceil(params.erp_vertical_offset_fraction * erp_wlen)), 2);
        otherwise
            error('Undefined vertical offset alignment method');
    end
    stacked_targets(:, tgt_indx, 1 : size(segment, 2)) = segment - mn;
end

erps_mean = zeros(num_ch, erp_wlen);
erps_var = zeros(num_ch, erp_wlen);
erps_rwa = zeros(num_ch, erp_wlen);
erps_rwmd = zeros(num_ch, erp_wlen);
for ch = 1 : num_ch
    block = squeeze(stacked_targets(ch, :, :));
    [erps_rwa(ch, :), ~, erps_rwmd(ch, :)] = robust_weighted_average(block);
    erps_mean(ch, :) = mean(block, 1);
    erps_var(ch, :) = var(block, [], 1);
end

if params.plot_results
    t = (0:sig_len-1)/fs;
    for i = 1:num_ch
        if mod(i, params.plots_per_figure) == 1 || params.plots_per_figure == 1
            figure;
        end
        subplot(params.plots_per_figure, 1, mod(i-1, params.plots_per_figure) + 1);
        plot(t, data(i,:),'b');
        hold on
        plot(t, data_detrended(i,:),'m');
        plot(t, baseline(i,:),'g');
        plot(t, data_den(i,:),'r');
        ylabel(num2str(i));
        grid;
        if mod(i, params.plots_per_figure) == 1 || params.plots_per_figure == 1
            title('Noisy vs. Denoised');
        end
        if mod(i,params.plots_per_figure) == 0 || params.plots_per_figure == 1
            xlabel('time(s)');
        end
    end

    tt = (0:erp_wlen-1)/fs;
    for i = 1:size(data, 1)
        figure
        hold on
        plot(tt, squeeze(stacked_targets(i, :, :))', 'color', [1, 0.6, 0.6]);
        h1 = plot(tt, erps_mean(i, :), 'k:', 'linewidth', 2);
        h2 = plot(tt, erps_rwa(i, :), 'k--', 'linewidth', 2);
        h3 = plot(tt, erps_rwmd(i, :), 'k', 'linewidth', 2);
        grid
        legend([h1, h2, h3], 'Mean', 'RWA', 'RWMD');
        set(gca, 'fontsize', 16)
        title(['Channel ', num2str(i)]);
    end

end

