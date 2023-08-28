function qtInt = qtIntGausFit_v1(ecg, fs, varargin)
% Estimate the Q-T interval by fitting Gaussian functions.
%
% Syntax:
% qtInt = qtIntGausFit_v1(ecg, fs)
% qtInt = qtIntGausFit_v1(ecg, fs, rpeaks)
% qtInt = qtIntGausFit_v1(ecg, fs, rpeaks, betaQ, betaT)
%
% Inputs:
%   ecg: Vector of single ECG data.
%   fs: Sampling frequency.
%   rpeaks (or hr): If vector, it is an index vector of R peak positions with the
%       same length as the ecg columns. rpeaks has 1 or true values at R peaks and 0
%       or false values elsewhere. If scalar, it is the average heart rate in Hz,
%       used in R peak detection.
%   betaQ & betaT: (Scalar) constant coefficient in (wave width) = 2 * beta * (Gaussian width), and 
%       wave offset/onset = Gaussian center +/- beta * (Gaussian width).
%
%   Note: Each of rpeaks, betaQ, and betaT can be left empty ([]), so the
%       default values will be used.
%
% Output:
%   qtInt: Q-T interval length.
%
% Reference:
%   Fattahi, Davood, and Reza Sameni. "Cramér-Rao Lower Bounds of
%   Model-Based Electrocardiogram Parameter Estimation." IEEE Transactions
%   on Signal Processing 70 (2022): 3181-3192.
%
% Revision History:
%   2021: First release
%
% Davood Fattahi (fattahi.d@gmail.com), 2021
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Initialize values
options = struct;
rPeaks = [];
beta.q = [];
beta.t = [];

%% Factor of onset/offset detection
if nargin >= 3
    rPeaks = varargin{1};
end
if nargin >= 4
    beta.q = varargin{2};
end
if nargin == 5
    beta.t = varargin{3};
end

% Segments of interest
soi.q(1) = -.080;
soi.q(2) = -.020;
soi.t(1) = .1;
soi.t(2) = .5;

% Initial points for Gaussian fitting optimization
p0 = [];

% Lower and upper bounds for Gaussian fitting optimization
lb = [];
ub = [];

% Optimization options
options = struct('SpecifyObjectiveGradient', true);

% R peaks
if size(rPeaks) == size(ecg)
    rPeaks = find(rPeaks);
end

% Estimate Q-T interval using Gaussian fitting
[~, ~, ~, ~, qtInt] = qtParamsGausFit_cl(ecg, fs, rPeaks, p0, soi, lb, ub, beta, options);
