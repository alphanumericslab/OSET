function [y, h] = quant_noise_wiener_filter(x, fs, nvar_est_mode, params)

if nargin < 3
    nvar_est_mode = 'min-spectral-power';
end

default_params = DefaultParameters();
if nargin < 4
    params = default_params;
end

if ~isfield(params, 'filter_len')
    params.filter_len = default_params.filter_len;
end
if ~isfield(params, 'innovation_filter_type')
    params.innovation_filter_type = default_params.innovation_filter_type;
end


switch nvar_est_mode
    case 'min-fft-power'
        nfft = min(1000, length(x));
        P = abs(fft(x, nfft)).^2 / nfft;
        lower_prctile = 5.0; % lower percentile threshold
        nvar = mean(P(P <= prctile(P, lower_prctile)));
    case 'min-spectral-power'
        nfft = min(1000, length(x));
        P = pwelch(x, ones(1, nfft)/nfft, round(0.75*nfft), nfft, fs , 'twosided'); % pwelch calculates one-sided spectra by default
        lower_prctile = 5.0; % lower percentile threshold
        nvar = mean(P(P <= prctile(P, lower_prctile)));
    case 'quantization-level' % to be used only on raw data collected from ADC without any filtering operations
        quanta = diff(sort(x)); % find all amplitude jumps
        quanta = quanta(quanta > 0); % remove the equal ones
        delta = mode(quanta); % ADC quantization level is the most prevalent
        nvar = delta.^2 / 12; % variance of uniformly distributed quantization noise
    case 'fix'
        nvar = params.nvar;
end

% Check and correct filter length for minimum-phase mode
if mod(params.filter_len, 2) == 0 && isequal(params.innovation_filter_type, 'MIN_PHASE')
    params.filter_len = params.filter_len + 1;
    warning(['Filter length needs to be odd in minimum-phase mode. New filter size: ', num2str(params.filter_len)]);
end

nfft = min(1000, length(x));
Px = pwelch(x, ones(1, nfft)/nfft, round(0.75*nfft), nfft, fs , 'twosided')'; % pwelch calculates one-sided spectra by default
Ps = Px - nvar;
Ps(Ps < eps) = eps;
H = Ps ./ Px;

f = fs * (0:length(H)-1)/length(H);
figure
plot(f, 10*log10(Ps))
hold on
plot(f, 10*log10(H))
grid

h0 = real(ifft(sqrt(H), params.filter_len, 2));
h0_zero = fftshift(h0, 2);
h0_zero = (h0_zero + h0_zero(:, end:-1:1)) / 2;

switch params.innovation_filter_type
    case 'LINEAR_PHASE'
        h = h0_zero;
    case 'MIN_PHASE'
        for ch = 1:N_channels
            r = conv(h0_zero(ch, :), h0_zero(ch, end:-1:1));
            h = firminphase(r);
        end
    otherwise
        error('Undefined innovations filter type');
end

% y = filter(h0_zero, 1, x);
y = conv(x, h0_zero);
delay = ceil((length(h0_zero) - 1) / 2); % Calculate the group delay
y = y(delay + 1 : delay + length(x)); % Compensate for the delay and keep the same length as x


end
% Function to generate default parameters
function params = DefaultParameters()
params.filter_len = 512;
params.innovation_filter_type = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE'
end
