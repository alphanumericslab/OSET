function [dipole, teta] = ecg_dipole_gen_state_space(N, fs, f, alphai, bi, tetai, teta0)
% 
% ecg_dipole_gen_state_space - Synthetic cardiac dipole generator using a state-space form of a sum of Gaussian model.
%   Refer to references for further details. Performs identical to dipole_gen_direct_sum.
%
% Usage:
%   [dipole, teta] = ecg_dipole_gen_state_space(N, fs, f, alphai, bi, tetai, teta0);
% 
% Inputs:
%   N: Signal length (in samples)
%   fs: Sampling rate
%   f: Heart rate (Hz)
%   alphai: Structure containing the amplitudes of Gaussian functions used
%       for modeling the x, y, and z coordinates of the cardiac dipole
%   bi: Structure containing the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   tetai: Structure containing the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   teta0: Initial phase of the synthetic dipole
%
% Outputs:
%   dipole: Structure containing the x, y, and z coordinates of the cardiac dipole
%   teta: Vector containing the dipole phase
%
% References:
%     - Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%       Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%       ECG Signals. In EURASIP Journal on Advances in Signal Processing (Vol.
%       2007, Issue 1). Springer Science and Business Media LLC.
%       https://doi.org/10.1155/2007/43407
%     - McSharry, P. E., Clifford, G. D., Tarassenko, L., & Smith, L. A.
%       (2003). A dynamical model for generating synthetic electrocardiogram
%       signals. In IEEE Transactions on Biomedical Engineering (Vol. 50, Issue
%       3, pp. 289–294). Institute of Electrical and Electronics Engineers
%       (IEEE). https://doi.org/10.1109/tbme.2003.808805
% 
% Revision History:
%   2006: First release
%   2023: Renamed from the deprecated version DipoleGenerator
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Constants
w = 2 * pi * f;
dt = 1 / fs;

teta = zeros(1, N);
X = zeros(1, N);
Y = zeros(1, N);
Z = zeros(1, N);

teta(1) = teta0;

for i = 1 : N-1
    % Update dipole phase
    teta(i+1) = mod(teta(i) + w * dt + pi, 2 * pi) - pi;

    % Calculate phase differences
    dtetaix = mod(teta(i) - tetai.x + pi, 2 * pi) - pi;
    dtetaiy = mod(teta(i) - tetai.y + pi, 2 * pi) - pi;
    dtetaiz = mod(teta(i) - tetai.z + pi, 2 * pi) - pi;

    if i == 1
        % Initial values for x, y, and z coordinates
        X(i) = sum(alphai.x .* exp(-dtetaix .^2 ./ (2 * bi.x .^ 2)));
        Y(i) = sum(alphai.y .* exp(-dtetaiy .^2 ./ (2 * bi.y .^ 2)));
        Z(i) = sum(alphai.z .* exp(-dtetaiz .^2 ./ (2 * bi.z .^ 2)));
    end

    % Update state variables
    X(i+1) = X(i) - dt * sum(w * alphai.x ./ (bi.x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2 * bi.x .^ 2)));
    Y(i+1) = Y(i) - dt * sum(w * alphai.y ./ (bi.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2 * bi.y .^ 2)));
    Z(i+1) = Z(i) - dt * sum(w * alphai.z ./ (bi.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2 * bi.z .^ 2)));
end

% Store x, y, and z coordinates in the dipole structure
dipole.x = X;
dipole.y = Y;
dipole.z = Z;
