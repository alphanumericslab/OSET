function [y, W, A] = sca_nsca_component_analysis(x, f1, f2, I, order)
% sca_nsca_component_analysis - Sort components by spectrally desired components vs nonstationary components.
%
% Usage:
%   [y, W, A] = sca_nsca_component_analysis(x, f1, f2, I, order)
%
% Inputs:
%   x: Input data array (channels x samples).
%   f1, f2: Normalized frequencies defining the frequency band of interest.
%   I: Indexes defining the desired/undesired components.
%   order: Sorting order of the spectral and nonstationary components ('forward' or 'backward').
%
% Outputs:
%   y: Sorted components based on spectrally desired components vs nonstationary components.
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
%   2023: Renamed from deprecated version SCA_NSCA.
%
% Reza Sameni, 2009-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Filter the signal and compute the covariance over the frequency band of interest
xx = bp_filter_dft(x, f1, f2);
B = xx * xx' / size(xx, 2);
B = (B + B') / 2;

% Compute the covariance over the indexes of interest
% C = cov(x(:, I)');
C = x(:, I) * x(:, I)' / length(I);
C = (C + C') / 2;

if isequal(order, 'forward')
    % Perform eigenvalue decomposition using Cholesky decomposition with B and C
    [V, D] = eig(B, C, 'chol');
elseif isequal(order, 'backward')
    % Perform eigenvalue decomposition using Cholesky decomposition with C and B
    [V, D] = eig(C, B, 'chol');
end

[~, II] = sort(diag(D), 'descend');

% Compute the extraction matrix
W = V(:, II)';
% Compute the mixing matrix (pseudo-inverse of the extraction matrix)
A = pinv(W);

% Extract the components by applying the extraction matrix to the input data
y = real(W * x);
end
