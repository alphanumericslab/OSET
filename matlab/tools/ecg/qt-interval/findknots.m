function k = findknots(x, nknots, varargin)
%
% Description:
%   This function finds suitable knots along the ECG signal for spline or
%   polynomial fitting. Starting with the first and last samples as initial
%   knots, it iteratively selects points with the maximum distance from the
%   reference line (the line passing through the knots) to set as new knots.
% 
% Usage:
%   k = findknots(x, nknots)
%   k = findknots(x, nknots, dist)
%
% Inputs:
%   x: A segment of the ECG signal.
%   nknots: The number of desired knots.
%   dist: Optional distance measure; accepts 'Euclidean' or 'Value'.
%
% Output:
%   k: Vector containing the found knot positions.
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%   on Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
%   2021: First release
%   2022:
%       1) Optional input 'nknots' added.
%       2) Unnecessary input parameters (depth and threshold) removed.
%       3) Bug of adjacent knots (one sample distance) resolved.
%       4) Line slope formula modified.
%       5) Optional inputs of value-distance and Euclidean distance added.
%
% Davood Fattahi (fattahi.d@gmail.com), 2021-2022
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Handle optional distance measure
if isempty(varargin)
    Ecd = true; Vald = false;
elseif strcmp(varargin{1}, 'Euclidean')
    Ecd = true; Vald = false;
elseif strcmp(varargin{1}, 'Value')
    Ecd = false; Vald = true;
end

% Check for valid number of knots
if nknots < 3
    error('Number of knots should be greater than 2.');
end

% Initialize knots matrix
x = x(:)';
k = nan(nknots - 1, 2);
k(1, 1) = 1;
k(1, 2) = size(x(:), 1);

% Iteratively find knots
for j = 2 : nknots - 1
    M = nan(j - 1, 1);
    I = nan(j - 1, 1);
    tt = cell(1, j - 1);
    
    % Calculate distances and positions
    for i = 1 : j - 1
        tt{i} = k(i, 1) : k(i, 2);
        s = x(tt{i});
        
        % Calculate distance based on the chosen method
        if Ecd
            m = (s(end) - s(1)) / (tt{i}(end) - tt{i}(1));
            [M(i), I(i)] = max(abs((-m * tt{i} + s) + m * tt{i}(1) - s(1)) / sqrt(m^2 + 1));
        elseif Vald
            line = s(1) + (0 : length(s) - 1) .* ((s(end) - s(1)) / (length(s) - 1));
            [M(i), I(i)] = max(abs(s - line)); % Value distance
        end
    end
    
    % Update knot positions
    [~, II] = max(M);
    k(j, 1) = tt{II}(I(II));
    k(j, 2) = k(II, 2);
    k(II, 2) = tt{II}(I(II));
end

% Ensure unique knot positions
k = unique(k(:));
