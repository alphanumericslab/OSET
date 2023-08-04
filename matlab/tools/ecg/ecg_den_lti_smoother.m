function y = ecg_den_lti_smoother(x, varargin)
% ecg_den_lti_smoother - Linear time invariant (LTI) ECG denoising based
%   on penalized least squares with smoothness priors (Tikhonov regularization).
%
% Usage:
%   y = ecg_den_lti_smoother(x, filter_param, SmoothnessFactor, mode, stable_pole_mag_th)
%
% Inputs:
%   x: Matrix of input ECG data (channels x samples)
%   filter_param: Either the filter impulse response (h) or the diff_order
%       parameter (optional). If diff_order is provided, the impulse
%       response h is calculated as a finite difference filter. Default is
%       a second order difference filter.
%   SmoothnessFactor: A non-negative scalar value representing the
%       smoothness penalty. Default is 1.0. It is used to control the
%       smoothness regularization level in the denoising process. Change
%       SmoothnessFactor by orders of magniture (in positive and negatize
%       powers of 10), to see the change.
%   mode: An integer representing the mode of operation (default: mode = 0):
%       - mode = 0 : y^* = argmin_y(|x - y| + lambda*|D*y + b|), fixed smoothness penalty (lambda)
%       - mode = 1: y^* = argmin_y(|D * y + b| + gamma*|x - y|), fixed MSE penalty (gamma)
%   stable_pole_mag_th: threshold of filter poles to be considered as stable (inside the unit circle). Default = 0.9999
%
% Note: Segment lengths should not be shorter than the filter length.
%
% Reference: See Section 4 in
%   R. Sameni, Online filtering using piecewise smoothness priors: Application
%   to normal and abnormal electrocardiogram denoising, Signal Processing,
%   Volume 133, 2017, Pages 52-63, https://doi.org/10.1016/j.sigpro.2016.10.019.
%
% Revision History:
%   2015: First release
%   2020: Modified to work with multichannel signals
%   2023: Renamed from deprecated version ECGSmoothnessPriorsDenoiserLTI
%
% Reza Sameni, 2015-2023
% The Open-Source Electrophysiological Toolbox

% use the given filter or calculate the filter impulse response
if nargin > 1 && ~isempty(varargin{1})
    filter_param = varargin{1};
    if length(filter_param) == 1
        diff_order = filter_param;
        h = diff([zeros(1, diff_order) 1 zeros(1, diff_order)], diff_order);
    else
        h = filter_param(:)';
    end
else
    diff_order = 2;
    h = diff([zeros(1, diff_order) 1 zeros(1, diff_order)], diff_order);
end

% SmoothnessFactor: default is 1
if nargin > 2 && ~isempty(varargin{2})
    SmoothnessFactor = varargin{2};
else
    SmoothnessFactor = 1;
    warning('Using default SmoothnessFactor = 1.0. Consider specifying the smoothness factor');
end

% mode: default is 0
if nargin > 3 && ~isempty(varargin{3})
    mode = varargin{3};
else
    mode = 0;
    warning('Using default mode = 0. Consider specifying the mode. See the function''s help for details.');
end

% mode: default is 0
if nargin > 4 && ~isempty(varargin{4})
    stable_pole_mag_th = varargin{4};
else
    stable_pole_mag_th = 0.9999;
end


phi = conv(h, h(end:-1:1));
midindex = (length(phi)+1)/2;
% See Section 4 of cited paper for the mathematical derivations
if mode == 0 % lambda-based
    gg = phi * SmoothnessFactor;
    gg(midindex) = gg(midindex) + 1;
elseif mode == 1 % gamma-based
    gg = phi;
    gg(midindex) = gg(midindex) + SmoothnessFactor;
end

r = roots(gg);
r_abs = abs(r);
I_causal = r_abs < stable_pole_mag_th;
r_causal = r(I_causal);
p_causal = poly(r_causal);
y = filtfilt(sum(p_causal), p_causal, x')'; % Note: the sqrt(SmoothnessFactor) which should be theoretically multiplied in x in this line is automatically compensated by the sum(p_causal) term
