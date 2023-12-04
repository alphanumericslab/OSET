function [Q, P] = PolarDecomposition(A, varargin)
% PolarDecomposition has been deprecated. Use polar_decomposition instead.
warning('PolarDecomposition has been deprecated. Use polar_decomposition instead.');
[Q, P] = polar_decomposition(A, varargin{:});