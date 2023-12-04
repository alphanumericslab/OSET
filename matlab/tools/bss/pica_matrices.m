function [A, B] = pica_matrices(x, peaks)
% pica_matrices - Calculates covariance and lagged-covariance matrices
%   required for the pseudo-periodic component analysis algorithm (PiCA).
%
% Usage:
%   [A, B] = pica_matrices(x, peaks)
%
% Inputs:
%   x: Matrix of input signals (rows represent channels, and columns represent samples)
%   peaks: Vector containing the indices of detected peaks in the input signal
%
% Outputs:
%   A: Covariance matrix between the peaks and the lagged time points
%   B: Covariance matrix of the data
%
% Reference:
%   R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
%       electrocardiogram decomposition using periodic component analysis. IEEE
%       Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version PiCAMatrixes1
%
% Reza Sameni, 2008-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Calculate the phase and discard the second output (not used)
[phase, ~] = phase_calculator(peaks);

% PM time calculation for the peaks sequence
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

% Calculate the covariance and lagged-covariance matrices for the peaks sequence
A = x(:, T0) * x(:, T1)';
B = x(:, T0) * x(:, T0)';
A = (A + A') / 2; % Symmetrize the covariance matrix
B = (B + B') / 2; % Symmetrize the lagged-covariance matrix
