function [y, p] = PolyFit(x, fs, N, varargin)
% PolyFit has been deprecated. Use polynomial_fit instead.
warning('PolyFit has been deprecated. Use polynomial_fit instead.');
[y, p] = polynomial_fit(x, fs, N, varargin{:});