function [p, varargout] = GausFit(t, x, p0, lb, ub, varargin)
%
% Fits a sum of Gaussian functions to the data.
%
% Usage:
% [p, varargout] = GausFit(t, x, p0, lb, ub)
% [p, varargout] = GausFit(t, x, p0, lb, ub, options)
% [p, varargout] = GausFit(t, x, p0, lb, ub, PrMu, PrCov, PstVar)
% [p, varargout] = GausFit(t, x, p0, lb, ub, PrMu, PrCov, PstVar, options)
%
% Inputs:
%   t: Time samples
%   x: Input data
%   p0: Initial parameter values. The first, second, and third rows respectively contain
%       Gaussian amplitudes (a), widths (b), and centers (c). Each column contains parameters
%       for a single Gaussian in the order [ai ; bi ; ci]. It can also be in vector form, like [a1 b1 c1 a2 b2 c2 ...]'.
%   lb: Lower bound for parameters
%   ub: Upper bound for parameters
%   options: Options structure for optimization problem (fields are the same as MATLAB's 'optimoptions')
%   PrMu: Mean of prior distribution
%   PrCov: Covariance matrix of prior information
%   PstVar: Variance of the additive noise (scalar)
%
% Outputs:
%   p: Estimated Gaussian parameters, with the same format as p0
%   varargout: Additional outputs depending on the optimization method
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%   on Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
% 2021: First release
%
% Davood Fattahi (fattahi.d@gmail.com), 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

p0 = p0(:);
t = t(:);
x = x(:);
lb = lb(:);
ub = ub(:);
NumGaus = length(p0) / 3;

% Validate the number of parameters
if round(NumGaus) ~= NumGaus
    error('Wrong number of parameters in p0!')
end

% Identify optimization method based on input arguments
if length(varargin) <= 1 % Maximum Likelihood
    ML = true;
    Bys = false;
elseif length(varargin) >= 3 % Bayesian
    ML = false;
    Bys = true;
    PrMu = varargin{1};
    PrMu = PrMu(:);
    PrCov = varargin{2};
    PstVar = varargin{3};
end

% Handle options input
if nargin == 5 || nargin == 8
    options0 = struct;
elseif nargin == 6 || nargin == 9
    options0 = varargin{end};
end

% Set default options if not provided
if isempty(options0)
    options0 = struct;
end

%%% Gaussian fit
if ML % Maximum Likelihood
    fun = @(p)(MlFun(t, p, x));
elseif Bys % Bayesian
    fun = @(p)(BysFun(t, p, x, PstVar, PrCov, PrMu));
end

% Set optimization options
options = optimoptions(@lsqnonlin, 'Display', 'off');
options = optionsmerge(options, options0);

% Perform optimization
[p, varargout{1:6}] = lsqnonlin(fun, p0, lb, ub, options);

% Sort Gaussian parameters based on ascending centers
p = reshape(p, 3, []);
[~, I] = sort(p(3, :));
p = p(:, I);
p = p(:);
end

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%%% Local functions

function [F, J] = MlFun(t, p, x)
    F = GausVal(t, p) - x;
    if nargout > 1
        J = GausGrad(t, p)';
    end
end

function [F, J] = BysFun(t, p, x, PstVar, PrCov, PrMu)
    F = ([(GausVal(t, p) - x) ./ sqrt(PstVar) ; sqrtm(PrCov) \ (p - PrMu)]);
    if nargout > 1
        J = [GausGrad(t, p)' ./ sqrt(PstVar) ; inv(sqrtm(PrCov))];
    end
end

function options = optionsmerge(options1, options2)
%%% Merge options2 into options1.
FN = fieldnames(options2);
for i = 1 : length(FN)
    eval(['options1.' FN{i} ' = options2.' FN{i} ';'])
end
options = options1;
end
