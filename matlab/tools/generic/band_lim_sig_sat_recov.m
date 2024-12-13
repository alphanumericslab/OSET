function y = band_lim_sig_sat_recov(x, th_bottom, th_top, n_itr, ff_low, ff_high, verbose)
% BAND_LIM_SIG_SAT_RECOV - Recovers a saturated signal by iterative band-limited filtering
%
% This function attempts to recover a saturated signal by applying
% iterative band-pass filtering while preserving the power of unsaturated
% parts of the signal. The algorithm roots in nonuniform sampling and has
% proven convergence properties for bandlimited signals with sparse
% saturations
%
% INPUTS:
%   x          : Input signal (vector)
%   th_bottom  : Lower saturation threshold
%   th_top     : Upper saturation threshold
%   n_itr      : Number of iterations for the recovery process
%   ff_low     : Low cutoff frequency for high-pass filtering normalized by
%                the sampling frequency (set to [] if no high-pass filter is desired)
%   ff_high    : High cutoff frequency for low-pass filtering normalized by
%                the sampling frequency
%
% OUTPUT:
%   y : Recovered signal
%
% Further Reading:
%   1- Marvasti, F. (1988). A unified approach to zero-crossings and
%       nonuniform sampling: Of single and multidimensional signal and systems.
%       Nonuniform. ISBN-10: 0961816708.
%   2- Marvasti, F. (Ed.). (2001). Nonuniform sampling: Theory and practice.
%       Kluwer Academic/Plenum Publishers. ISBN-10: 0306464454
%
% Reza Sameni, 2024
% The Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET

if nargin < 7
    verbose = false;
end

% Check if the signal is saturated
if max(x) <= th_top && min(x) >= th_bottom
    warning('Signal is not saturated; returning input signal unchanged');
    y = x;
    return
end

% Identify the indexes of the unsaturated parts of the signal
maintained_indexes = find(x <= th_top & x >= th_bottom);

% Compute the power of the unsaturated part of the signal
pre_filter_unsaturated_signal_power = mean(x(maintained_indexes).^2);

% Initialize
y = x;

% Iterate
for k = 1 : n_itr
    % Preserve the unsaturated sampled of the signal in each iteration
    y(maintained_indexes) = x(maintained_indexes);

    % Apply a low-pass filter with zero phase
    y = lp_filter_zero_phase(y, ff_high);

    % Apply a high-pass filter
    if ~isempty(ff_low)
        y = y - lp_filter_zero_phase(y, ff_low);
    end

    % Normalize the power to the unsaturated portions' of the original signal
    y = y / mean(y(maintained_indexes).^2) * pre_filter_unsaturated_signal_power;

    if verbose
        disp(['Iteration = ', num2str(k)]);
    end
end
end
