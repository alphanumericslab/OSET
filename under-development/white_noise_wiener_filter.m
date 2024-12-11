function [y, h] = white_noise_wiener_filter(x, fs, nvar_est_mode, params)
% WHITE_NOISE_WIENER_FILTER - Designs and applies a Wiener filter to suppress
%   white noise in a signal with estimated or provided noise variance
%
% Inputs:
%   x              - Input signal
%   fs             - Sampling frequency (Hz)
%   nvar_est_mode  - Noise variance estimation mode:
%                    'min-spectral-power' (default), 'min-fft-power',
%                    'quantization-level', 'fix'
%   params         - Structure containing optional parameters for filtering.
%                    See DefaultParameters() for default values.
%
% Outputs:
%   y              - Filtered signal
%   h              - Impulse response of the Wiener filter
%
% Reference: Papoulis, A., & Pillai, S. U. (2002). Probability, random
%   variables,  and stochastic processes (4th ed.). McGraw-Hill.
% 
% Reza Sameni, 2024
% The Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET

% Check the number of input arguments and set default noise variance estimation mode if necessary
if nargin < 3
    nvar_est_mode = 'min-spectral-power';
end

% Load default parameters if not provided
default_params = DefaultParameters();

if nargin < 4
    params = default_params;
end

% Ensure all fields are present in the params structure
fields = {'filtering_method', 'plot_response', 'filter_len', 'lower_prctile', ...
    'max_nfft', 'spectral_est_overlap', 'innovation_filter_type', ...
    'output_normalization_method'};
for i = 1:length(fields)
    if ~isfield(params, fields{i})
        params.(fields{i}) = default_params.(fields{i});
    end
end

% Estimate noise variance based on the specified mode
switch nvar_est_mode
    case 'min-fft-power'
        % Use FFT power spectrum for noise variance estimation
        nfft = min(params.max_nfft, length(x));
        P = abs(fft(x, nfft)).^2 / nfft;
        nvar = mean(P(P <= prctile(P, params.lower_prctile)));
    case 'min-spectral-power'
        % Use Welch's method for noise variance estimation
        nfft = min(params.max_nfft, length(x));
        P = pwelch(x, ones(1, nfft)/nfft, floor(params.spectral_est_overlap*nfft), nfft, fs, 'twosided');
        nvar = mean(P(P <= prctile(P, params.lower_prctile)));
    case 'quantization-level'
        % Estimate noise from input quantization levels (only applicable to unprocessed ADC/quantized data)
        quanta = diff(sort(x)); % Find amplitude jumps
        quanta = quanta(quanta > 0); % Remove zeros
        delta = mode(quanta); % ADC quantization level
        nvar = delta.^2 / 12; % Variance of uniform quantization noise (is a standard formula). See slide 348 onward here: https://www.researchgate.net/publication/325581786_Digital_Systems_Design_Course_Slides
    case 'fix'
        % Use a fixed noise variance provided as input in params
        nvar = params.nvar;
end

% Adjust filter length for minimum-phase mode
if mod(params.filter_len, 2) == 0 && isequal(params.innovation_filter_type, 'MIN_PHASE')
    params.filter_len = params.filter_len + 1;
    warning(['Filter length adjusted to odd value for minimum-phase mode: ', num2str(params.filter_len)]);
end

% Design the Wiener filter
nfft = min(params.max_nfft, length(x));
Px = pwelch(x, hamming(nfft), round(params.spectral_est_overlap*nfft), nfft, fs, 'twosided')'; % Power spectrum
Ps = Px - nvar; % Signal power spectrum
Ps(Ps < eps) = eps; % Ensure non-negative power
H = Ps ./ Px; % Wiener filter frequency response

% Construct filter impulse response
h0 = real(ifft(sqrt(H), params.filter_len, 2));
h0_zero = fftshift(h0, 2); % Zero-phase version
h0_zero = h0_zero + h0_zero(:, end:-1:1);

% Apply the specified innovation filter type
switch params.innovation_filter_type
    case 'LINEAR_PHASE'
        h = h0_zero;
    case 'MIN_PHASE' % minimum-phase filter design from impulse response
        for ch = 1:size(h0_zero, 1)
            r = conv(h0_zero(ch, :), h0_zero(ch, end:-1:1));
            h = firminphase(r); % Minimum-phase design
        end
    otherwise
        error('Undefined innovation filter type');
end

% Apply the Wiener filter to the input signal
switch params.filtering_method
    case 'filter'
        y = filter(h, 1, x);
    case 'conv'
        y = conv(x, h);
        delay = ceil((length(h) - 1) / 2); % Group delay
        y = y(delay + 1 : delay + length(x)); % Compensate for delay
    otherwise
        error('Undefined filtering method');
end

% Normalize the output signal if required
switch params.output_normalization_method
    case 'none'
        % No normalization
    case 'match-input-power'
        y = y / sqrt(sum(abs(h).^2)); % Match input signal power
    case 'match-dc-gain'
        y = y / sqrt(sum(h)); % Match DC gain
end

% Plot filter response if enabled
if params.plot_response
    f = fs * (0:length(H)-1)/length(H); % Frequency vector
    figure;
    hold on;
    plot(f, 10*log10(Px), 'DisplayName', 'Px (input spectra)');
    plot(f, 10*log10(Ps), 'DisplayName', 'Ps (signal spectra)');
    plot(f, 10*log10(nvar*ones(1, length(H))), 'DisplayName', 'Pn (noise spectra)');
    xlim([f(1), f(end)])
    legend show;
    title('Input, output and noise spectra')
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('Wiener Filter Response');
    set(gca, 'fontsize', 16)

    figure
    stem(h)
    title("Wiener filter's impulse response")
    xlabel('n');
    grid
    set(gca, 'fontsize', 16)
end

end

% Function to generate default parameters
function params = DefaultParameters()
% DefaultParameters - Returns a structure with default parameter values.
params.lower_prctile = 1.0;                  % Lower percentile threshold
params.max_nfft = 1024;                      % Maximum FFT size
params.spectral_est_overlap = 0.75;         % Overlap for spectral estimation
params.filter_len = 512;                    % Wiener filter length
params.innovation_filter_type = 'LINEAR_PHASE'; % 'LINEAR_PHASE' or 'MIN_PHASE'
params.filtering_method = 'conv';           % Filtering method ('conv' or 'filter')
params.output_normalization_method = 'none'; % Output normalization method
params.plot_response = true;                % Plot filter response flag
end
