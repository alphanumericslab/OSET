function [Q, P] = PolarDecomposition(A, varargin)
%
% Calculates the polar decomposition of a square (complex) matrix A:
%    A = QP in the right mode
%    A = PQ in the left mode
%  where Q is a unitary matrix and P is a positive-semidefinite Hermitian
%  matrix.
%
% Usage:
% [Q, P] = PolarDecomposition(A, mode);
% 
% Inputs:
% A: A square (complex) matrix
% mode: 'right' (default) or 'left'
%
% Outputs:
% Q: a unitary matrix
% P: a positive-semidefinite Hermitian matrix
%
% Reza Sameni, Jun 2019
%
% Open Source Electrophysiological Toolbox, version 3.14, June 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% reza.sameni@gmail.com
%

[m, n] = size(A);
if(m ~= n)
    error('Input matrix should be square');
end

if(nargin == 1)
    side = 'right';
elseif(nargin == 2)
    side = varargin{1};
    if(~isequal(side, 'right') && ~isequal(side, 'left'))
        error('Unknown mode');
    end
else
    error('Up to two inputs are acceptable');
end

[U,S,V] = svd(A);
Q = U*V';
if(isequal(side, 'right'))
    P = V*S*V';
else
    P = U*S*U';
end
