function [Q, P] = polar_decomposition(A, varargin)
%
% Calculates the polar decomposition of a square (complex) matrix A:
%       A = QP in the right mode
%       A = PQ in the left mode
%   where Q is a unitary matrix and P is a positive-semidefinite Hermitian
%   matrix.
%
% Usage:
%   [Q, P] = polar_decomposition(A, mode);
%
% Inputs:
%   A: A square (complex) matrix
%   mode: 'right' (default) or 'left'
%
% Outputs:
%   Q: a unitary matrix
%   P: a positive-semidefinite Hermitian matrix
%
% Revision History:
%   2019: First release
%   2023: Renamed from deprecated version PolarDecomposition.
%
% Reza Sameni, 2019-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check if the input matrix is square
[m, n] = size(A);
if m ~= n
    error('Input matrix should be square');
end

% Determine the mode (right or left) based on input arguments
if nargin > 1 && ~isempty(varargin{1})
    side = varargin{1};
    if ~isequal(side, 'right') && ~isequal(side, 'left')
        error('Unknown mode');
    end
else
    side = 'right'; % Default to 'right' mode
end

% Perform the singular value decomposition
[U, S, V] = svd(A);

% Calculate the unitary matrix Q
Q = U * V';

% Calculate the positive-semidefinite Hermitian matrix P based on the mode
switch side
    case 'right'
        P = V * S * V';
    case 'left'
        P = U * S * U';
end
