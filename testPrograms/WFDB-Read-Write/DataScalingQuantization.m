function [data_quant, bias, gains, max_amp, snr, snr_med, num_quant_levs] = DataScalingQuantization(data, method, params)
%
% DataScalingQuantization - Perform data quantization and analysis on multi-channel data.
%
% Syntax:
%   [data_quant, bias, gains, max_amp, snr, snr_med, num_quant_levs] = DataScalingQuantization(data, method, params)
%
% Input Arguments:
%   - data: A matrix of size (num_ch, num_samples) representing the input data, where num_ch is the number of channels and num_samples is the number of samples.
%   - method: A string specifying the data scaling/saturation method. Possible values are:
%             - 'max_scale': No saturation or clipping.
%             - 'tanh_sat': Saturate the channels at k_sigma_sat times the channel STD using the hyperbolic tangent function.
%             - 'clip': Clip the channels at a specified clip_level.
%   - params: A structure containing additional parameters for the quantization process:
%             - params.quant_bits: Number of quantization bits (1 to 64).
%             - params.k_sigma_sat: Scaling factor for saturation (used in 'tanh_sat' method).
%             - params.clip_level: Clip level (used in 'clip' method).
%             - params.plot_mode: Plotting mode for visualization. Possible values are:
%               - 'samples': Plot original data, saturated/clipped data, and reconstructed data.
%               - 'errors': Plot the difference between original and saturated/clipped data, and the difference between original and quantized data.
%               - 'noplots': Does not plot
%             - params.remove_mean: remove the mean of each channel before
%               digitization and report the bias (true or false)
%
% Output Arguments:
%   - data_quant: A matrix of the same size as the input data, representing the quantized version of the data.
%   - bias: A matrix of the same size as the input data, storing the biases (channel means) that were removed if params.remove_mean == true.
%   - gains: A vector containing the maximum gain values for each channel.
%   - max_amp: Maximum quantized output amplitude: 2^(params.quant_bits-1) - 1
%   - snr: A vector containing the Signal-to-Noise Ratio (SNR) in decibels for each channel.
%   - snr_med: A vector containing the SNR in decibels for each channel, using the median error power.
%   - num_quant_levs: A vector containing the number of quantization levels for each channel in the original data.
%
% Example:
%   data = randn(4, 1000); % Generate random multi-channel data
%   method = 'tanh_sat';
%   params.quant_bits = 8;
%   params.k_sigma_sat = 3;
%   params.plot_mode = 'samples';
%   [data_quant, bias, gains, max_amp, snr, snr_med, num_quant_levs] = DataScalingQuantization(data, method, params);
%
% Note:
%   1- An approximation of the original double-precision samples from the quantized version can be calculated as follows (per channel ch):
%   data_approx(ch, :) = bias(ch) + double(data_quant(ch, :)) / gains(ch)
%
%   2- params.quant_bits can be integers between 1 to 64. The quantized output will be type-casted upwards to: 8, 16, 32, or 64
%
% Reza Sameni, May 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

num_ch = size(data, 1); % number of channels

if params.quant_bits < 1 || params.quant_bits > 64 || mod(params.quant_bits,1) ~= 0
    error(['Unsupported number of quantization bits: ', num2str(params.quant_bits)]);
elseif params.quant_bits >= 1 && params.quant_bits <= 8
    quant_dtype = 'int8';
elseif params.quant_bits <= 16
    quant_dtype = 'int16';
elseif params.quant_bits <= 32
    quant_dtype = 'int32';
else
    quant_dtype = 'int64';
end

% Find the number of quantization levels of the original data
num_quant_levs = zeros(1, size(data, 1));
for ch = 1 : num_ch
    num_quant_levs(ch) = length(find(diff(sort(data(ch, :))))) + 1;
end

% Remove channel means and save the biases
bias = zeros(1, size(data, 1));
if params.remove_mean % remove mean if true
    for ch = 1 : num_ch
        bias(ch) = mean(data(ch, :));
        data(ch, :) = data(ch, :) - bias(ch);
    end
end

% Saturate/clip the amplitudes
data_sat = zeros(size(data));
switch method
    case 'max_scale' % No saturation/clipping
        for ch = 1 : num_ch
            data_sat(ch, :) = data(ch, :);
        end
    case 'tanh_sat' % Saturate the channels at k_sigma_sat times the channel STD
        for ch = 1 : num_ch
            gamma = params.k_sigma_sat * std(data(ch, :));
            if gamma > max(abs(data(ch, :))) % do not saturate if channel max is already below k_sigma_sat times the STD
                data_sat(ch, :) = data(ch, :);
            else
                data_sat(ch, :) = gamma * tanh(data(ch, :) / gamma);
            end
        end
    case 'clip' % Saturate the channels at k_sigma_sat times the channel STD
        for ch = 1 : num_ch
            data_sat(ch, :) = data(ch, :);
            I_clip = data_sat(ch, :) > params.clip_level;
            data_sat(ch, I_clip) = params.clip_level;
        end
    otherwise
        error('Undefined data scaling/saturation method.');
end

gains = zeros(1, num_ch);
max_amp = 2^(params.quant_bits-1) - 1; % maximum values after quantization. Note: WFDB considers -2^(quant_bits-1) as NAN
data_quant = zeros(size(data), quant_dtype); % the quantized version of the data
for ch = 1 : num_ch
    gains(ch) = max_amp / max(abs(data_sat(ch, :)));
    data_quant(ch, :) = cast(round(gains(ch) * data(ch, :)), quant_dtype);
end

snr = zeros(1, num_ch);
snr_med = zeros(1, num_ch);
data_approx = zeros(size(data));
for ch = 1 : num_ch
    data_approx(ch, :) = double(data_quant(ch, :)) / gains(ch);
    er = data(ch, :) - data_approx(ch, :);
    snr(ch) = 10*log10(mean(data(ch, :).^2) / mean(er.^2));
    snr_med(ch) = 10*log10(mean(data(ch, :).^2) / median(er.^2));
end

switch params.plot_mode
    case 'samples'
        for ch = 1 : num_ch
            lgnd = {};
            figure
            hold on
            plot(data(ch, :)); lgnd = cat(1, lgnd, 'data');
            plot(data_sat(ch, :)); lgnd = cat(1, lgnd, 'saturated/clipped');
            plot(data_approx(ch, :)); lgnd = cat(1, lgnd, 'reconstructed');
            title(['Ch#: ', num2str(ch)]);
            legend(lgnd, 'interpreter', 'none');
            ylabel('Amplitude');
            xlabel('Samples');
            set(gca, 'fontsize', 18);
            grid
        end
    case 'errors'
        for ch = 1 : num_ch
            lgnd = {};
            figure
            hold on
            plot(data(ch, :) - data_sat(ch, :)); lgnd = cat(1, lgnd, 'data-saturated/clipped');
            plot(data(ch, :) - data_approx(ch, :)); lgnd = cat(1, lgnd, 'data - quantized');
            title(['Ch#: ', num2str(ch)]);
            legend(lgnd, 'interpreter', 'none');
            ylabel('Amplitude');
            xlabel('Samples');
            set(gca, 'fontsize', 18);
            grid
        end
    case 'noplots'
        % Don't plot
    otherwise
        error('Undefined plot mode');
end
