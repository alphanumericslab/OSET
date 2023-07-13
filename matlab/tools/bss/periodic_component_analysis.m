function [y, W, A] = periodic_component_analysis(x, peaks1, varargin)
% periodic_component_analysis - Pseudo-Periodic Component Analysis (PiCA).
%   Extracts and ranks sources in order of resemblance with
%   pseudo-periodic events (such as ECG) with one or two sets of cyclic
%   peaks, from multichannel data
%
% Usage:
%   [y, W, A] = periodic_component_analysis(x, peaks1, peaks2, mode)
%
% Inputs:
%   x: Input data array (channels x samples).
%   peaks1: Vector of first signal peaks impulse train (R-wave locations
%       for ECG signals).
%   peaks2 (optional): Vector of second signal peaks impulse train (R-wave
%       locations for ECG signals).
%   mode (optional): Mode indicating the behavior of the algorithm:
%        0: Single-peak case, components ranked based on
%           resemblance to peaks1 periodicity. (default mode if only one peak is selected)
%        1: Two-peak case, components ranked based on resemblance to
%           peaks1 - peaks2. (default mode if two peaks are given as input)
%        2: Two-peak case, blending from similarity with peaks1 to peaks2.
%
% Outputs:
%   y: Pseudo-periodic components ranked in order of resemblance with the
%       first to second (if available) desired signals.
%   W: Decomposition (separation) matrix.
%   A: Mixing matrix.
%
% Note: If the signal peaks2 is not available, the components are ranked
%   from the most similar to the least similar to the signal corresponding to
%   peaks1.
%
% Reference:
%   R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
%       electrocardiogram decomposition using periodic component analysis. IEEE
%       Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2008: First release 2023: Renamed from deprecated version PiCA
%
% Reza Sameni, 2008-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin > 2 && ~isempty(varargin{1})
    peaks2 = varargin{1};
    flag = 1; % two-peak case
else
    flag = 0; % single-peak case
end

if nargin > 3 % two-peak case in blending mode
    flag = varargin{2};
end

% Calculate PM time
[T0, T1] = synchronous_phase_samples(peaks1);

% Calculate covariance or lagged-covariance matrices
A = x(:, T0) * x(:, T1)';
B = x(:, T0) * x(:, T0)';

% Make matrices symmetric, if otherwise
A = (A + A') / 2;
B = (B + B') / 2;

if flag == 0 % Single peaks case
    [V, D] = eig(A, B, 'chol');
elseif flag == 1 % Two peak sets case
    [T0, T1] = synchronous_phase_samples(peaks2);
    AA = x(:, T0) * x(:, T1)';
    AA = (AA + AA') / 2;
    [V, D] = eig(A - AA, B, 'chol');
elseif flag == 2 % Two peak sets case blending from similarity with the first to the second set
    [T0, T1] = synchronous_phase_samples(peaks2);
    AA = x(:, T0) * x(:, T1)';
    AA = (AA + AA') / 2;
    [V, D] = eig(A, AA, 'chol');
end

[~, I] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order

W = V(:, I)'; % Extraction matrix
A = pinv(W); % Mixing matrix (pseudo-inverse of the extraction matrix)

y = real(W * x); % Extract spectral components by applying the extraction matrix to the input data