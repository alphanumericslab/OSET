function [y,W,A] = spectral_component_analysis_dft(x,fl,fu)
% spectral_component_analysis_dft - Spectral component analysis (SCA) using
% a Fourier domain approach. Extracts and ranks linear mixtures of
% multichannel data with maximal energy in a given frequency band.
%
% Usage:
%   [y,W,A] = spectral_component_analysis_dft(x,fl,fu)
%
% Inputs:
%   x: Input data array (channels x samples).
%   fl: lower cutoff frequency normalized by sampling frequency
%   fu: upper cutoff frequency normalized by sampling frequency
%
% Outputs:
%   y: Extracted spectral components ranked by their energy in the
%       frequency band of interest
%   W: Extraction (separation) matrix
%   A: Mixing matrix
%
% Reference: - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2009). A
%   deflation procedure for subspace decomposition. IEEE Transactions on
%   Signal Processing, 58(4), 2363-2374.
%
% Revision History:
%   2009: First release
%   2023: Renamed from deprecated version SCA. corrected a frequency
%       binning bug
%
% Reza Sameni, 2009-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

L2 = size(x,2);

freq_indexes = max(round(L2*fl), 1) : min(round(L2*fu), L2); % indexes of the frequency range of interest (half band only)

X = fft(x, L2, 2); % Fourier transform of the input data
Bf = X(:, freq_indexes) * X(:, freq_indexes)'; % Compute the spectral covariance matrix
B = 2 * real(Bf); % Multiply by 2 to compensate for counting only half the signal power in the frequency range of interest
B = (B + B') / 2; % Make the matrix symmetric

C = x * x'; % Compute the covariance matrix
C = (C + C') / 2; % Make the matrix symmetric

[V, D] = eig(B, C, 'chol'); % Perform eigenvalue decomposition using Cholesky decomposition

[~, II] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order

W = V(:, II)'; % Extraction matrix
A = pinv(W); % Mixing matrix (pseudo-inverse of the extraction matrix)

y = real(W * x); % Extract spectral components by applying the extraction matrix to the input data
end
