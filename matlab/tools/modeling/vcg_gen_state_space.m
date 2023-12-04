function [vcg, phi] = vcg_gen_state_space(N, fs, f, alpha, b, theta, theta0)
% 
% vcg_gen_state_space - Synthetic cardiac dipole/vectorcardiogram (VCG) generator using a state-space form of a sum of Gaussian model.
%   Refer to references for further details. Performs identical to dipole_gen_direct_sum.
%
% Usage:
%   [vcg, phi] = vcg_gen_state_space(N, fs, f, alpha, b, theta, theta0);
% 
% Inputs:
%   N: Signal length (in samples)
%   fs: Sampling rate
%   f: Heart rate (Hz)
%   alpha: Structure containing the amplitudes of Gaussian functions used
%       for modeling the x, y, and z coordinates of the cardiac dipole/VCG
%   b: Structure containing the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole/VCG
%   theta: Structure containing the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole/VCG
%   theta0: Initial phase of the synthetic dipole/VCG
%
% Outputs:
%   vcg: Structure containing the x, y, and z coordinates of the cardiac dipole/VCG
%   phi: Vector containing the dipole/VCG phase
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

phi = zeros(1, N);
X = zeros(1, N);
Y = zeros(1, N);
Z = zeros(1, N);

phi(1) = theta0;

for i = 1 : N-1
    % Update dipole/VCG phase
    phi(i+1) = mod(phi(i) + w * dt + pi, 2 * pi) - pi;

    % Calculate phase differences
    dtetaix = mod(phi(i) - theta.x + pi, 2 * pi) - pi;
    dtetaiy = mod(phi(i) - theta.y + pi, 2 * pi) - pi;
    dtetaiz = mod(phi(i) - theta.z + pi, 2 * pi) - pi;

    if i == 1
        % Initial values for x, y, and z coordinates
        X(i) = sum(alpha.x .* exp(-dtetaix .^2 ./ (2 * b.x .^ 2)));
        Y(i) = sum(alpha.y .* exp(-dtetaiy .^2 ./ (2 * b.y .^ 2)));
        Z(i) = sum(alpha.z .* exp(-dtetaiz .^2 ./ (2 * b.z .^ 2)));
    end

    % Update state variables
    X(i+1) = X(i) - dt * sum(w * alpha.x ./ (b.x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2 * b.x .^ 2)));
    Y(i+1) = Y(i) - dt * sum(w * alpha.y ./ (b.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2 * b.y .^ 2)));
    Z(i+1) = Z(i) - dt * sum(w * alpha.z ./ (b.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2 * b.z .^ 2)));
end

% Store x, y, and z coordinates in the dipole/VCG structure
vcg.x = X;
vcg.y = Y;
vcg.z = Z;
