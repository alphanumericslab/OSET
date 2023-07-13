function [y, W, A] = spectral_component_analysis_bp(x, fl, fu, varargin)
% spectral_component_analysis_bp - Spectral component analysis (SCA) using
%   time-domain bandpass filtering approach. Extracts and ranks linear mixtures of
%   multichannel data with maximal energy in a given frequency band.
%
% Usage:
%   [y, W, A] = spectral_component_analysis_bp(x, fl, fu, method)
%
% Inputs:
%   x: Input data array (channels x samples).
%   fl: Lower cutoff frequency normalized by the sampling frequency.
%   fu: Upper cutoff frequency normalized by the sampling frequency.
%   method (optional): Filtering method to be used ('DFT_FILTER' or 'DESIGN_FILTER').
%
% Outputs:
%   y: Extracted spectral components ranked by their energy in the frequency band of interest.
%   W: Extraction (separation) matrix.
%   A: Mixing matrix.
% 
% Note: Depending on the frequency specifications the filter design may be
%   unstable when method = 'DESIGN_FILTER'
%
% Reference: Sameni, R., Jutten, C., & Shamsollahi, M. B. (2009). A
%   deflation procedure for subspace decomposition. IEEE Transactions on
%   Signal Processing, 58(4), 2363-2374.
%
% Revision History:
%   2009: First release
%   2023: Renamed from deprecated version SCAmodified. Corrected a frequency binning bug.
%
% Reza Sameni, 2009-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if nargin > 3 && ~isempty(varargin{1})
    method = varargin{1};
else
    method = 'DFT_FILTER';
end

switch method
    case 'DFT_FILTER'
        x_filtered = bp_filter_dft(x, fl, fu); % Filter the input data using a DFT-based approach
    case 'DESIGN_FILTER'
        h = fdesign.bandpass('N,Fc1,Fc2,Ast1,Ap,Ast2', 80, fl, fu, 50, 0.5, 50, 1); % Design a bandpass filter
        Hd = design(h);
        hdf1 = convert(Hd,'df1'); % Convert the filter to direct form I
        x_filtered = filter(hdf1, x, 2); % Apply the designed filter to the input data
    otherwise
        error('Undefined filtering method');
end

B = x_filtered * x_filtered'; % Compute the covariance matrix of the filtered data
C = x * x'; % Compute the covariance matrix of the original data

C = (C + C') / 2; % Make the covariance matrix symmetric
B = (B + B') / 2; % Make the covariance matrix of the filtered data symmetric

[V,D] = eig(B, C, 'chol'); % Perform eigenvalue decomposition using Cholesky decomposition

[~,II] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order

W = V(:, II)'; % Extraction matrix
A = pinv(W); % Mixing matrix (pseudo-inverse of the extraction matrix)

y = real(W * x); % Extract spectral components by applying the extraction matrix to the input data