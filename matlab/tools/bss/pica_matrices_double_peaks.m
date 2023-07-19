function [A, AA, B] = pica_matrices_double_peaks(x, peaks, peaks2)
% pica_matrices_double_peaks - Calculates covariance and lagged-covariance matrices
%   required for the pseudo-periodic component analysis algorithm (PiCA). 
%
% Usage:
%   [A, AA, B] = pica_matrices_double_peaks(x, peaks, peaks2)
%
% Inputs:
%   x: Matrix of input signals (rows represent channels, and columns represent samples)
%   peaks: Vector containing the indices of detected peaks in the input signal (for the first sequence)
%   peaks2: Vector containing the indices of detected peaks in the input signal (for the second sequence)
%
% Outputs:
%   A: Covariance matrix between the first sequence (peaks) and the lagged time points
%   AA: Covariance matrix between the second sequence (peaks2) and the second sequence's time points
%   B: Covariance matrix of the data
%
% Reference:
%   R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
%       electrocardiogram decomposition using periodic component analysis. IEEE
%       Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version PiCAMatrixes
%
% Reza Sameni, 2008-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Calculate the ECG phase
[phase, ~] = phase_calculator(peaks);
[phase2, ~] = phase_calculator(peaks2);

% PM time calculation for the first peaks sequence
J = find(peaks);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1) - min(n1);

N0 = length(x);
T1 = zeros(length(peaks) - prd - wlen, 1);
NN = length(T1);
for t = 1:NN
    df = abs(phase(t) - phase(max(t + prd - wlen, 1):min(t + prd + wlen, N0)));
    [~, I] = min(df);
    T1(t) = t + prd + I - wlen - 1;
end
T1 = max(T1, 1);
T1 = min(T1, N0);
T0 = 1:NN;

% Calculate the covariance and lagged-covariance matrices for the first peaks sequence
A = x(:, T0) * x(:, T1)';
B = x(:, T0) * x(:, T0)';
A = (A + A') / 2; % Symmetrize the covariance matrix
B = (B + B') / 2; % Symmetrize the lagged-covariance matrix

% PM time calculation for the second peaks sequence
J = find(peaks2);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1) - min(n1);

T1 = zeros(length(peaks2) - prd - wlen, 1);
NN = length(T1);
for t = 1:NN
    df = abs(phase2(t) - phase2(max(t + prd - wlen, 1):min(t + prd + wlen, N0)));
    [~, I] = min(df);
    T1(t) = t + prd + I - wlen - 1;
end
T1 = max(T1, 1);
T1 = min(T1, N0);
T0 = 1:NN;

% Calculate the covariance matrix for the second peaks sequence
AA = x(:, T0) * x(:, T1)';
AA = (AA + AA') / 2; % Symmetrize the covariance matrix
