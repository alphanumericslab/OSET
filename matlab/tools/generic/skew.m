function [skw, m, sd] = skew(data)
% skew - Skewness of a matrix data (over the second dimension).
%
% Syntax: [skw, m, sd] = skew(data)
%
% Inputs:
%   data: N channels x T samples data matrix.
%
% Outputs:
%   skw: Skewness vector (N x 1).
%   m: Mean vector (N x 1).
%   sd: Standard deviation vector (N x 1).
%
%
%   Revision History:
%       2019: First release
%       2023: Added documentation
%
%   Reza Sameni, 2019-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

m = mean(data, 2);      % Compute the mean along the second dimension
sd = std(data, [], 2);  % Compute the standard deviation along the second dimension
m3 = mean(data.^3, 2);  % Compute the mean of data^3 along the second dimension

skw = (m3 - 3*m.*sd.^2 - m.^3) ./ sd.^3;  % Calculate the skewness
