function [y, p] = polynomial_fit(x, fs, N, varargin)
% Polynomial fit over a signal segment using least square or
%   pseudo-inversion
%
% [y, p] = polynomial_fit(x, fs, N, method)
%
% Inputs:
%   x: the (noisy) input signal
%   fs: sampling frequency
%   N: polynomial order (>= 1)
%   method: optional argument to specify the method ('LS': least squares or 'PINV': pseudo inverse) for
%       solving the linear equation. Default = 'LS'
%
% Outputs:
%   y: fitted mode
%   p: polynomial coefficients
%
%   Revision History:
%       2011: First release
%       2023: Documented and extended to LS and renamed from deprecated version PolyFit()
%
%   Reza Sameni, 2011-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

x = x(:)'; % Ensure x is a row vector

if N < 0
    error('polynomial order must be >= 0');
end

if nargin > 3
    if isequal(varargin{1}, 'LS') || isequal(varargin{1}, 'PINV')
        method = varargin{1};
    else
        error('Unknown method');
    end
else
    method = 'LS';
end

L = length(x);
t = (0 : L-1) / fs;

T = zeros(N + 1, L);
TX = zeros(N + 1, 1);
for i = 1 : N + 1
    T(i, :) = t.^(i-1); % Construct the matrix T with powers of t
    TX(i) = x * T(i,:)'; % Calculate the inner product of x and the i-th row of T
end

switch method
    case 'LS'
        p = (T * T') \ TX; % Solve the linear equation T * T' * p = TX for p using least squares
    case 'PINV'
        p = pinv(T * T') * TX; % Solve the linear equation T * T' * p = TX for p using the pseudoinverse
end

y = p' * T; % Calculate the fitted mode y using the polynomial coefficients p and T
