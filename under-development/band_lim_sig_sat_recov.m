function y = band_lim_sig_sat_recov(x, th_bottom, th_top, n_itr, ff_low, ff_high, th_bottom_end, th_top_end, verbose)
% BAND_LIM_SIG_SAT_RECOV - Recovers a saturated signal by iterative band-limited filtering
%
% This function attempts to recover a saturated signal by applying
% iterative band-pass filtering while preserving the power of unsaturated
% parts of the signal. The algorithm roots in nonuniform sampling theory 
% for bandlimited signals with sparse saturations.
%
% INPUTS:
%   x              : Input signal (vector)
%   th_bottom      : Lower saturation threshold
%   th_top         : Upper saturation threshold
%   n_itr          : Number of iterations for the recovery process
%   ff_low         : Low cutoff frequency for high-pass filtering, normalized
%                    by the sampling frequency (set to [] if no high-pass
%                    filter is desired)
%   ff_high        : High cutoff frequency for low-pass filtering, normalized
%                    by the sampling frequency
%   th_bottom_end  : Lower bound for saturated values after recovery.
%                    If not specified, it is set to 
%                    `th_bottom - 0.25 * abs(th_bottom)` by default.
%   th_top_end     : Upper bound for saturated values after recovery.
%                    If not specified, it is set to 
%                    `th_top + 0.25 * abs(th_top)` by default.
%   verbose        : Boolean flag indicating whether to display progress
%                    messages during iterations. Default is `false`.
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

% The default factor beyond max/min that the signal should never exceed
default_upper_lower_extreme_scale_factor = 0.25;

if nargin < 7 || isempty(th_bottom_end)
    th_bottom_end = th_bottom - default_upper_lower_extreme_scale_factor*abs(th_bottom);
end

if nargin < 8 || isempty(th_top_end)
    th_top_end = th_top + default_upper_lower_extreme_scale_factor*abs(th_top);
end

if nargin < 9
    verbose = false;
end

% Check if the signal is saturated
if max(x) <= th_top && min(x) >= th_bottom
    warning('Signal is not saturated; returning input signal unchanged');
    y = x;
    return
end

% Identify the indexes of the unsaturated parts of the signal
maintained_indexes = x <= th_top & x >= th_bottom;
% saturated_indexes = ~ maintained_indexes;

% Compute the power of the unsaturated part of the signal
pre_filter_unsaturated_signal_power = std(x(maintained_indexes));
% pre_filter_unsaturated_signal_power = std(x);

% Remove mean
x_dc = mean(x);
x = x - x_dc;

% Initialize
y = x;

% set the satrurated segments to the extreme max/min range
y(y > th_top) = th_top_end;
y(y < th_bottom) = th_bottom_end;

% Iterate
for k = 1 : n_itr
    % Preserve the unsaturated sampled of the signal in each iteration
    y(maintained_indexes) = x(maintained_indexes);

    % Make sure the output does not exceed the defined extreme cases
    y(y > th_top_end) = th_top_end;
    y(y < th_bottom_end) = th_bottom_end;

    % Apply a low-pass filter with zero phase
    y_bp = lp_filter_zero_phase(y, ff_high);

    % Apply a high-pass filter
    if ~isempty(ff_low)
        y_bp = y_bp - lp_filter_zero_phase(y_bp, ff_low);
    end

    % Normalize the power to maintain the power of the unsaturated portions' of the original signal
    y = y_bp / std(y_bp(maintained_indexes)) * pre_filter_unsaturated_signal_power;
    % Normalize to the input power
    % y = y / std(y) * pre_filter_unsaturated_signal_power;

    if verbose
        disp(['Iteration = ', num2str(k)]);
    end
end

% Add back the signal mean
y = y + x_dc;

end
