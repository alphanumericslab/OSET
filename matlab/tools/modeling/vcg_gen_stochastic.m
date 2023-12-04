function [vcg, phi]= vcg_gen_stochastic(N, fs, f, f_deviations, alpha, delta_alpha, b, delta_b, theta, delta_theta, teta0)
% vcg_gen_stochastic - Synthetic cardiac dipole/vectorcardiogram (VCG) generator with beat-wise stochastic deviations
%
% [vcg, phi]= vcg_gen_stochastic(N, fs, f, f_deviations, alpha, delta_alpha, b, delta_b, theta, delta_theta, teta0)
%
% This function generates a synthetic cardiac dipole/VCG signal with stochastic
%   deviations at each beat. The deviations affect heart rate, dipole/VCG amplitudes,
%   Gaussian width, and Gaussian phase.
%
% Inputs:
%   N: signal length
%   fs: sampling rate
%   f: average heart rate (Hz)
%   f_deviations: percentage of beat-wise heart rate deviations (Hz)
%   alpha: structure containing the amplitudes of Gaussian functions for
%       modeling x, y, and z coordinates of the cardiac dipole/VCG
%   delta_alpha: percentage of amplitude deviations added per beat
%   b: structure containing the widths of Gaussian functions for modeling
%       x, y, and z coordinates of the cardiac dipole/VCG
%   delta_b: percentage of Gaussian width deviations added per beat
%   theta: structure containing the phase of Gaussian functions for modeling
%       x, y, and z coordinates of the cardiac dipole/VCG
%   delta_theta: percentage of Gaussian center deviations added per beat
%   teta0: initial phase of the synthetic dipole/VCG
%
% Outputs:
%   vcg: structure containing the x, y, and z coordinates of the cardiac dipole/VCG
%   phi: vector containing the dipole/VCG phase
%
% References:
%   - Clifford, G. D., Nemati, S., & Sameni, R. (2010). An artificial vector
%     model for generating abnormal electrocardiographic rhythms. Physiological
%     Measurement, 31(5), 595–609.
%   - Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%     Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%     ECG Signals. EURASIP Journal on Advances in Signal Processing, 2007(1).
%   - McSharry, P. E., Clifford, G. D., Tarassenko, L., & Smith, L. A.
%     (2003). A dynamical model for generating synthetic electrocardiogram signals.
%     IEEE Transactions on Biomedical Engineering, 50(3), 289–294.
%
% Revision History:
%   2022: First release
%   2023: Renamed from the deprecated version DipoleGeneratorStochastic
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

w = 2*pi*f; % Angular frequency
dt = 1/fs; % Time step

phi = zeros(1,N); % Initialize dipole/VCG phase vector
X = zeros(1,N); % Initialize x-coordinate vector
Y = zeros(1,N); % Initialize y-coordinate vector
Z = zeros(1,N); % Initialize z-coordinate vector

phi(1) = teta0; % Set initial dipole/VCG phase
d_alphai.x = alpha.x; % Initialize dipole/VCG amplitude deviations
d_alphai.y = alpha.y;
d_alphai.z = alpha.z;

d_tetai.x = theta.x; % Initialize dipole/VCG phase deviations
d_tetai.y = theta.y;
d_tetai.z = theta.z;

d_bi.x = b.x; % Initialize Gaussian width deviations
d_bi.y = b.y;
d_bi.z = b.z;

for i = 1 : N-1
    dtetaix = mod(phi(i) - d_tetai.x + pi , 2*pi) - pi; % Phase difference for x-coordinate
    dtetaiy = mod(phi(i) - d_tetai.y + pi , 2*pi) - pi; % Phase difference for y-coordinate
    dtetaiz = mod(phi(i) - d_tetai.z + pi , 2*pi) - pi; % Phase difference for z-coordinate

    if i == 1
        % Compute initial values for the first iteration
        X(i) = sum(d_alphai.x .* exp(-dtetaix .^2 ./ (2*d_bi.x .^ 2)));
        Y(i) = sum(d_alphai.y .* exp(-dtetaiy .^2 ./ (2*d_bi.y .^ 2)));
        Z(i) = sum(d_alphai.z .* exp(-dtetaiz .^2 ./ (2*d_bi.z .^ 2)));
    end

    % Update x, y, and z state variables
    X(i+1) = X(i) - dt*sum(w * d_alphai.x ./ (d_bi.x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* d_bi.x .^ 2)));
    Y(i+1) = Y(i) - dt*sum(w * d_alphai.y ./ (d_bi.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* d_bi.y .^ 2)));
    Z(i+1) = Z(i) - dt*sum(w * d_alphai.z ./ (d_bi.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* d_bi.z .^ 2)));

    % Update parameters for the next beat
    phi(i+1) = phi(i) + w*dt;
    if phi(i+1) > pi % Beat transition
        phi(i+1) = phi(i+1) - 2*pi;
        
        % Apply stochastic deviations to parameters
        d_alphai.x = alpha.x .* (1 + (rand(1, length(alpha.x)) - 0.5) * delta_alpha);
        d_alphai.y = alpha.y .* (1 + (rand(1, length(alpha.y)) - 0.5) * delta_alpha);
        d_alphai.z = alpha.z .* (1 + (rand(1, length(alpha.z)) - 0.5) * delta_alpha);

        d_tetai.x = theta.x .* (1 + (rand(1, length(theta.x)) - 0.5) * delta_theta);
        d_tetai.y = theta.y .* (1 + (rand(1, length(theta.y)) - 0.5) * delta_theta);
        d_tetai.z = theta.z .* (1 + (rand(1, length(theta.z)) - 0.5) * delta_theta);

        d_bi.x = b.x .* max(0, (1 + (rand(1, length(b.x)) - 0.5) * delta_b));
        d_bi.y = b.y .* max(0, (1 + (rand(1, length(b.y)) - 0.5) * delta_b));
        d_bi.z = b.z .* max(0, (1 + (rand(1, length(b.z)) - 0.5) * delta_b));
        
        % Update angular frequency with stochastic deviation
        w = 2 * pi * f * max(0, (1 + (rand - 0.5) * f_deviations));
    end
end

% Store the dipole/VCG coordinates in the output structure
vcg.x = X;
vcg.y = Y;
vcg.z = Z;
