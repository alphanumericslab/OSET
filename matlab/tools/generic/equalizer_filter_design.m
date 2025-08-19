function [h, H2] = equalizer_filter_design(Sx, Sy, params)
%EQUALIZER_FILTER_DESIGN Design a frequency-domain equalizer from power spectra.
%
%   Syntax:
%       [h, H2] = equalizer_filter_design(Sx, Sy)
%       [h, H2] = equalizer_filter_design(Sx, Sy, params)
%
%   Description:
%       This function designs an equalizer filter from the given source and
%       target power spectra. It optionally smooths the spectra, computes the
%       equalizer magnitude-squared response, and then derives a time-domain
%       filter with either linear-phase or minimum-phase characteristics.
%       An optional plot of the spectra and filter impulse response can be
%       generated.
%
%   Inputs:
%       Sx      - Source power spectrum (vector, length N)
%       Sy      - Target power spectrum (vector, length N)
%       params  - (Optional) Structure of parameters. If omitted or incomplete,
%                 defaults from DefaultParameters() are used. Fields include:
%
%         .fs                        - (double) Sampling frequency in Hz
%         .filter_len                - (integer) Length of time-domain filter
%         .lambda                    - (double) Regularization parameter for
%                                      Tikhonov smoothing
%         .smooth_spectrum           - (logical) Whether to smooth input spectra
%         .innovation_filter_type    - (string) Filter phase type
%                                      ('LINEAR_PHASE', 'MIN_PHASE')
%         .plot_results              - (logical) Whether to plot results
%
%   Outputs:
%       h   - Designed equalizer impulse response (vector)
%       H2  - Equalizer squared magnitude frequency response
%
%   Example:
%       % Example with defaults
%       [h, H2] = equalizer_filter_design(Sx, Sy);
%
%   See also: fft, ifft, fftshift, firminphase

    % ---- Default parameter assignment -------------------------------------
    def = DefaultParameters();
    if nargin < 3 || isempty(params)
        params = def;
    else
        fields_def = fieldnames(def);
        for k = 1:numel(fields_def)
            f = fields_def{k};
            if ~isfield(params, f) || isempty(params.(f))
                params.(f) = def.(f);
            end
        end
    end

    % ---- Input checks -----------------------------------------------------
    Sx = Sx(:)'; % row vector
    Sy = Sy(:)'; % row vector
    if length(Sx) ~= length(Sy)
        error('Sx and Sy must have the same length.');
    end
    spectral_len = length(Sx);

    % ---- Optional smoothing -----------------------------------------------
    if params.smooth_spectrum
        wlen = 7;
        threshold = 0.05;
        [Sx_smoothed, ~] = conditional_median_filter(Sx, wlen, threshold, 0.3);
        % Sx_smoothed = tikhonov_regularization(Sx, 2, params.lambda);
        Sx_smoothed(Sx_smoothed < params.epsilon) = params.epsilon;
        Sx_smoothed = (Sx_smoothed + flipud(Sx_smoothed)) / 2;

        wlen = 7;
        threshold = 0.05;
        [Sy_smoothed, ~] = conditional_median_filter(Sy, wlen, threshold, 0.3);
        % Sy_smoothed = tikhonov_regularization(Sy, 2, params.lambda);
        Sy_smoothed(Sy_smoothed < params.epsilon) = params.epsilon;
        Sy_smoothed = (Sy_smoothed + flipud(Sy_smoothed)) / 2;
    else
        Sx_smoothed = Sx;
        Sy_smoothed = Sy;
    end

    % ---- Equalizer magnitude-squared --------------------------------------
    H2 = Sy_smoothed ./ Sx_smoothed;

    % ---- Time-domain prototype filter -------------------------------------
    h0      = real(ifft(sqrt(abs(H2)), params.filter_len));
    h0_zero = fftshift(h0);
    h0_zero = (h0_zero + flipud(h0_zero)) / 2; % enforce even symmetry

    % ---- Phase type selection ---------------------------------------------
    switch upper(params.innovation_filter_type)
        case 'LINEAR_PHASE'
            h = h0_zero;

        case 'MIN_PHASE'
            r = conv(h0_zero, flipud(h0_zero));
            h = firminphase(r);

        otherwise
            error('Undefined innovation filter type: %s', params.innovation_filter_type);
    end

    % ---- Optional plotting ------------------------------------------------
    if params.plot_results
        ff = params.fs * (0:spectral_len-1) / spectral_len;

        figure
        subplot(1,2,1)
        lgnd = {};
        plot(ff - params.fs/2, fftshift(10*log10(abs(Sx))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sx';
        plot(ff - params.fs/2, fftshift(10*log10(abs(Sx_smoothed))), 'linewidth', 2); lgnd{end+1} = 'Sx_smoothed';
        plot(ff - params.fs/2, fftshift(10*log10(abs(Sy))), 'linewidth', 3); lgnd{end+1} = 'Sy';
        plot(ff - params.fs/2, fftshift(10*log10(abs(Sy_smoothed))), 'linewidth', 2); lgnd{end+1} = 'Sy_smoothed';
        plot(ff - params.fs/2, fftshift(10*log10(abs(H2))), 'linewidth', 3); lgnd{end+1} = 'H^2';
        xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
        legend(lgnd, 'interpreter', 'none')
        grid on
        title('Original (Sx) vs Target (Sy) spectra and Equalizer Squared Magnitude (H^2)');
        set(gca, 'fontsize', 18)

        subplot(1,2,2)
        stem(h, 'filled'); grid on
        xlabel('n'); ylabel('Impulse Response'); set(gca, 'fontsize', 18)

        figure
        freqz(h, 1, 1000, "whole", params.fs);
        set(gca, 'fontsize', 18)
    end
end

% ========================= Default Parameters =============================
function params = DefaultParameters()
    params.epsilon                   = eps;
    params.fs                        = 1.0;
    params.filter_len                = 512;
    params.lambda                     = 10000.0;
    params.smooth_spectrum           = true;
    params.innovation_filter_type    = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE'
    params.plot_results              = true;
end
