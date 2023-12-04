function baseline = baseline_filter(x_raw, baseline_removal_method, params)
% baseline_filter - Different approaches for baseline wander estimation using OSET functions.
%
% Syntax: baseline = baseline_filter(x_raw, baseline_removal_method, params)
%
% Inputs:
%   x_raw: Raw input signal (channels x time).
%   baseline_removal_method: Method for baseline removal.
%   params: Parameters required for the selected method.
%
% Output:
%   baseline: Estimated baseline wander signal.
%
%   Revision History:
%       2021: First release
%       2023: Renamed from deprecated version BaselineEstimator()
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

switch baseline_removal_method
    case 'BYPASS' % Bypass baseline wander (zero baseline).
        baseline = zeros(size(x_raw));
    case 'LP' % Low-pass filter.
        baseline = lp_filter_zero_phase(x_raw, params.fc / params.fs);
    case 'MNMN' % Two-stage moving average.
        baseline = baseline_sliding_window_twice(x_raw, params.wlen1, params.wlen2, 'mn');
    case 'MDMN' % Two-stage moving median plus moving average.
        % Apply a two-stage moving-median (md) and moving-average (mn) baseline detector.
        bl1 = baseline_sliding_window(x_raw, params.wlen1, 'md');
        baseline = baseline_sliding_window(bl1, params.wlen2, 'mn');
    case 'MDMDMDMN' % Three median filters followed by a moving average.
        bl1 = baseline_sliding_window(x_raw, params.wlen1, 'md');
        bl2 = baseline_sliding_window(x_raw, params.wlen2, 'md');
        bl3 = baseline_sliding_window(x_raw, params.wlen3, 'md');
        baseline = baseline_sliding_window((bl1 + bl2 + bl3) / 3, params.wlen4, 'mn');
    case 'TIKHONOV' % Tikhonov regularization.
        baseline = tikhonov_regularization(x_raw, params.DiffOrder, params.lambda);
    otherwise
        error('Unknown baseline removal method');
end
