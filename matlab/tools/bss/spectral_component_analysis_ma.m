function [y, W, A] = spectral_component_analysis_ma(x, f0, bw, order)
% spectral_component_analysis_ma - Spectral component analysis (SCA) using multi-stage moving average (MA) filtering approach. Extracts and ranks linear mixtures of multichannel data with maximal energy in a given frequency band.
%
% Usage:
%   [y, W, A] = spectral_component_analysis_ma(x, f0, bw, order)
%
% Inputs:
%   x: Input data array (channels x samples).
%   f0: Center frequency of the frequency band of interest, normalized by
%       the sampling frequency
%   bw: Bandwidth of the frequency band of interest, normalized by
%       the sampling frequency
%   order: Order of the multi-stage moving average filter.
%
% Outputs:
%   y: Extracted spectral components ranked by their energy in the frequency band of interest.
%   W: Extraction (separation) matrix.
%   A: Mixing matrix.
%
% Reference:
%   Sameni, R., Jutten, C., & Shamsollahi, M. B. (2009). A deflation
%   procedure for subspace decomposition. IEEE Transactions on Signal
%   Processing, 58(4), 2363-2374.
%
% Revision History:
%   2009: First release
%   2023: Renamed from deprecated version SCA2.
%
% Reza Sameni, 2009-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Perform multi-stage moving average (MA) filtering
x_filtered = bp_filter_complex_ma(x, f0, bw, order);

% Compute the covariance matrix of the filtered data
B = x_filtered * x_filtered';

% Compute the covariance matrix of the original data
C = x * x';

% Make both covariance matrices symmetric
C = (C + C') / 2;
B = (B + B') / 2;

% Perform eigenvalue decomposition using the Cholesky decomposition
[V, D] = eig(B, C, 'chol');

% Sort eigenvalues in descending order
[~, II] = sort(diag(D), 'descend');

% Compute the extraction matrix
W = V(:, II)';

% Compute the mixing matrix (pseudo-inverse of the extraction matrix)
A = pinv(W);

% Extract spectral components by applying the extraction matrix to the input data
y = real(W * x);
end
