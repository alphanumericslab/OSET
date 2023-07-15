function y = lp_filter_zero_phase(x, fc)
% lp_filter_zero_phase - Second-order zero-phase Lowpass filter.
%
% Syntax: y = lp_filter_zero_phase(x, fc)
%
% Inputs:
%   x: Vector or matrix of input data (channels x samples).
%   fc: -3dB cut-off frequency normalized by the sampling frequency.
%
% Output:
%   y: Vector or matrix of filtered data (channels x samples).
%
%   Revision History:
%       2006: First release
%       2023: Renamed from deprecated version LPFilter()
% 
%   Read more:  Mitra, S. (2010). Digital signal processing (4th ed.)
%               New York, NY: McGraw-Hill Professional.
%
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if fc >= 1
    error('fc should be smaller than 1');
end

k = 0.7071; % Cut-off value of 1/sqrt(2) or -6dB amplitude attenuation
alpha = (1 - k*cos(2*pi*fc) - sqrt(2*k*(1 - cos(2*pi*fc)) - k^2*sin(2*pi*fc)^2)) / (1 - k); % analytically derived. See ref. Mitra (2010)
y = filtfilt(1 - alpha, [1, -alpha], x')';
