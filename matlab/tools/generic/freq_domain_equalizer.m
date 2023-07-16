function x_filtered = freq_domain_equalizer(x, method, params, varargin)
% freq_domain_equalizer - Frequency domain equalization
% This function implements several frequency domain equalizers that alter
% the frequency-domain magnitude response while preserving the phase.
%
% Syntax:
%   x_filtered = freq_domain_equalizer(x, method, params, plot_flag)
%
% Inputs:
%   x: Input multichannel signal, where each row represents a channel and
%       each column represents a sample.
%   method: A string specifying the equalization method.
%       Available options:
%       - 'GLOBAL_EQ': Applies Tikhonov regularization for global
%                      frequency domain equalization.
%       - 'MED_PERCENTILE': Uses the median filter for smoothing based on
%                          percentile thresholds.
%       - 'MED_MSE': Uses the median filter for smoothing based on MSE.
%       - 'MEAN_MSE': Uses the mean filter (moving average) for smoothing based on MSE.
%   params: A structure containing the parameters for the selected method.
%       The required fields for each method are as follows:
%       - 'GLOBAL_EQ':
%         - 'diff_order_or_filter_coefs': Order of the Tikhonov regularization
%           model for frequency domain equalization.
%         - 'lambda': Regularization parameter for the Tikhonov regularization.
%       See help for: tikhonov_regularization
%       - 'MED_PERCENTILE':
%         - 'half_wlen': Half-width of the median filter window in frequency
%           samples (integer value). Larger values provide more smoothing
%           but may reduce signal resolution.
%         - 'filtering_perctile': Percentile (in the range 0 to 100) used to
%           determine the threshold for spectral point filtering. Spectral
%           amplitudes above this percentile are considered for outlier
%           replacement. Suggested range: 5 to 10 (percentiles above this threshold
%           are conditionally smoothed in the frequency domain).
%         - 'outlier_mse_percentage': Outlier neighborhood percentage used to
%           determine the outlier threshold. The magnitude difference between
%           the filtered and original Fourier transform is multiplied by
%           this percentage to determine outlier points. Suggested range: 5 to 20
%           (signal dependent).
%       - 'MED_MSE':
%         - 'half_wlen': Half-width of the median filter window in frequency
%           samples (integer value). Larger values provide more smoothing
%           but may reduce signal resolution.
%         - 'filtering_perctile': Percentile (in the range 0 to 100) used to
%           determine the threshold for spectral point filtering. Spectral
%           amplitudes above this percentile are considered for outlier
%           replacement. Suggested range: 5 to 10 (percentiles above this threshold
%           are conditionally smoothed in the frequency domain).
%         - 'outlier_mse_percentage': Outlier neighborhood percentage used to
%           determine the outlier threshold. The magnitude difference between
%           the filtered and original Fourier transform is multiplied by
%           this percentage to determine outlier points. Suggested range: 5 to 20
%           (signal dependent).
%       - 'MEAN_MSE':
%         - 'half_wlen': Half-width of the mean filter window in frequency
%           samples (integer value). Larger values provide more smoothing
%           but may reduce signal resolution.
%         - 'filtering_perctile': Percentile (in the range 0 to 100) used to
%           determine the threshold for spectral point filtering. Spectral
%           amplitudes above this percentile are considered for outlier
%           replacement. Suggested range: 5 to 10 (percentiles above this threshold
%           are conditionally smoothed in the frequency domain).
%         - 'outlier_mse_percentage': Outlier neighborhood percentage used to
%           determine the outlier threshold. The magnitude difference between
%           the filtered and original Fourier transform is multiplied by
%           this percentage to determine outlier points. Suggested range: 5 to 20
%           (signal dependent).
%   plot_flag (optional, default: 0): A flag (0 or 1) to indicate whether to plot the
%       results or not. Default is 0 (no plot).
%
% Outputs:
%   x_filtered: Filtered multichannel signal with frequency domain spikes
%   replaced by smoothed values based on the adaptive filtering approach.
%
% Description:
%   The function first applies a Fourier transform on the input
%   multichannel signal 'x'. It then applies one of the specified equalization
%   methods in the frequency domain with the specified parameters to smooth the
%   Fourier magnitude at each frequency point, while preserving the
%   phase. Finally, an inverse Fourier transform is applied to obtain the
%   filtered multichannel signal in the time domain.
%
% Example:
%   x = randn(8, 1000); % Example input multichannel noise signal
%   params.half_wlen = 5;
%   params.filtering_perctile = 5;
%   params.outlier_mse_percentage = 10;
%   plot_flag = 1;
%   x_filtered = freq_domain_equalizer(x, 'MED_MSE', params, plot_flag);
%
% Revision History:
%   2023: First release
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin > 3 && ~isempty(varargin{1})
    plot_flag = varargin{1};
else
    plot_flag = 0;
end

N = size(x, 1);
T = size(x, 2);
X = fft(x, T, 2); % Fourier transform

X_shifted = fftshift(X, 2); % swap the left and right parts of the spectra (DC is shifted to the center)
X_shifted_abs = abs(X_shifted); % the Fourier transform magnitude

switch method
    case 'GLOBAL_EQ'
        if isfield(params, 'diff_order_or_filter_coefs')
            diff_order_or_filter_coefs = params.diff_order_or_filter_coefs;
        else
            diff_order_or_filter_coefs = 2;
            warning('Undefined Tikhonov regularization model order set to 2');
        end
        if isfield(params, 'lambda')
            lambda = params.lambda;
        else
            error('params.lambda undefined');
        end
        X_shifted_abs_refined = tikhonov_regularization(X_shifted_abs, diff_order_or_filter_coefs, lambda);
    case 'MED_PERCENTILE'
        if isfield(params, 'filtering_perctile')
            filtering_perctile = params.filtering_perctile;
        else
            error('params.filtering_perctile undefined');
        end
        if isfield(params, 'half_wlen')
            half_wlen = params.half_wlen;
        else
            error('params.half_wlen undefined');
        end
        if isfield(params, 'outlier_mse_percentage')
            outlier_mse_percentage = params.outlier_mse_percentage;
        else
            error('params.outlier_mse_percentage undefined');
        end
        X_shifted_abs_filtered = zeros(size(X_shifted));
        I_filtered_indexes = false(size(X_shifted));
        filtering_threshold = prctile(X_shifted_abs, filtering_perctile, 2); % filters only amplitudes above this percentile
        low = min([outlier_mse_percentage, 100 - outlier_mse_percentage]);
        high = max([outlier_mse_percentage, 100 - outlier_mse_percentage]);
        % Median filter
        for t = 1 : T
            indexes = max(t - half_wlen, 1) : min(t + half_wlen, T);
            pctl = prctile(X_shifted_abs(:, indexes), [low, high], 2);
            I_filtered_indexes(:, t) = (X_shifted_abs(:, t) < pctl(:, 1) | X_shifted_abs(:, t) > pctl(:, 2)) ...
                & X_shifted_abs(:, t) > filtering_threshold;
            X_shifted_abs_filtered(:, t) = median(X_shifted_abs(:, indexes), 2); % the smoothed magnitude of the Fourier transform
        end

        % Conditional frequency domain smoothing
        X_shifted_abs_refined = X_shifted_abs;
        X_shifted_abs_refined(I_filtered_indexes) = X_shifted_abs_filtered(I_filtered_indexes); % smoothed spectra only replaces outlier spectral points (corresponding to the outliers)
    case 'MED_MSE'
        if isfield(params, 'filtering_perctile')
            filtering_perctile = params.filtering_perctile;
        else
            error('params.filtering_perctile undefined');
        end
        if isfield(params, 'half_wlen')
            half_wlen = params.half_wlen;
        else
            error('params.half_wlen undefined');
        end
        if isfield(params, 'outlier_mse_percentage')
            outlier_mse_percentage = params.outlier_mse_percentage;
        else
            error('params.outlier_mse_percentage undefined');
        end
        X_shifted_abs_filtered = zeros(size(X_shifted));
        % Median filter
        for t = 1 : T
            indexes = max(t - half_wlen, 1) : min(t + half_wlen, T);
            X_shifted_abs_filtered(:, t) = median(X_shifted_abs(:, indexes), 2); % the smoothed magnitude of the Fourier transform
        end

        % Conditional frequency domain smoothing
        original_smoothed_diff = X_shifted_abs - X_shifted_abs_filtered; % absolute error between filtered and original spectral magnitudes
        filtering_threshold = prctile(X_shifted_abs, filtering_perctile, 2); % filters only amplitudes above this percentile
        I_filtered_indexes = abs(original_smoothed_diff) > outlier_mse_percentage * X_shifted_abs_filtered & X_shifted_abs > repmat(filtering_threshold, 1, T); % spectral points corresponding to the grid

        X_shifted_abs_refined = X_shifted_abs;
        X_shifted_abs_refined(I_filtered_indexes) = X_shifted_abs_filtered(I_filtered_indexes); % smoothed spectra only replaces outlier spectral points (corresponding to the outliers)
    case 'MEAN_MSE'
        if isfield(params, 'filtering_perctile')
            filtering_perctile = params.filtering_perctile;
        else
            error('params.filtering_perctile undefined');
        end
        if isfield(params, 'half_wlen')
            half_wlen = params.half_wlen;
        else
            error('params.half_wlen undefined');
        end
        if isfield(params, 'outlier_mse_percentage')
            outlier_mse_percentage = params.outlier_mse_percentage;
        else
            error('params.outlier_mse_percentage undefined');
        end
        X_shifted_abs_filtered = zeros(size(X_shifted));
        % Median filter
        for t = 1 : T
            indexes = max(t - half_wlen, 1) : min(t + half_wlen, T);
            X_shifted_abs_filtered(:, t) = mean(X_shifted_abs(:, indexes), 2); % the smoothed magnitude of the Fourier transform
        end

        % Conditional frequency domain smoothing
        original_smoothed_diff = X_shifted_abs - X_shifted_abs_filtered; % absolute error between filtered and original spectral magnitudes
        filtering_threshold = prctile(X_shifted_abs, filtering_perctile, 2); % filters only amplitudes above this percentile
        I_filtered_indexes = abs(original_smoothed_diff) > outlier_mse_percentage * X_shifted_abs_filtered & X_shifted_abs > repmat(filtering_threshold, 1, T); % spectral points corresponding to the grid

        X_shifted_abs_refined = X_shifted_abs;
        X_shifted_abs_refined(I_filtered_indexes) = X_shifted_abs_filtered(I_filtered_indexes); % smoothed spectra only replaces outlier spectral points (corresponding to the outliers)
end

X_abs_refined = ifftshift(X_shifted_abs_refined, 2);
X_abs_refined(:, 2:end) = (X_abs_refined(:, 2:end) + X_abs_refined(:, end:-1:2)) / 2; % make Fourier magnitude even symmetrical, excluding the DC component (frequency index = 1)

X_shifted_phase = angle(X); % the Fourier transform phase
X_filtered =  X_abs_refined .* exp(1j*X_shifted_phase); % the magnitude smoothed Fourier transform (the phase remains unchanged)
x_filtered = real(ifft(X_filtered, T, 2)); % Inverse Fourier transform

% Plot results if in verbose mode
if plot_flag
    for ch = 1 : N
        figure
        plot(x(ch, :));
        hold on
        plot(x_filtered(ch, :));
        title(['Ch #', num2str(ch)]);
        legend('Original', 'Filtered');
        xlabel('Samples');
        ylabel('Amplitude');
        grid
    end

    for ch = 1 : N
        f = (0:T-1)/T;
        figure
        plot(f, 20*log10(abs(X(ch, :))));
        hold on
        plot(f, 20*log10(abs(X_filtered(ch, :))));
        legend('Original', 'Filtered');
        title(['Ch #', num2str(ch), ' spectrum']);
        xlabel('Normalized frequency');
        ylabel('Magnitude (dB)');
        grid
    end
end
