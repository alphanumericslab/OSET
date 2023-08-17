function [dipole, teta]= ecg_dipole_gen_stochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
% ecg_dipole_gen_stochastic - Synthetic cardiac dipole generator with beat-wise stochastic deviations
%
% [dipole teta]= ecg_dipole_gen_stochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
%
% This function generates a synthetic cardiac dipole signal with stochastic
%   deviations at each beat. The deviations affect heart rate, dipole amplitudes,
%   Gaussian width, and Gaussian phase.
%
% Inputs:
%   N: signal length
%   fs: sampling rate
%   f: average heart rate (Hz)
%   f_deviations: percentage of beat-wise heart rate deviations (Hz)
%   alphai: structure containing the amplitudes of Gaussian functions for
%       modeling x, y, and z coordinates of the cardiac dipole
%   delta_alphai: percentage of amplitude deviations added per beat
%   bi: structure containing the widths of Gaussian functions for modeling
%       x, y, and z coordinates of the cardiac dipole
%   delta_bi: percentage of Gaussian width deviations added per beat
%   tetai: structure containing the phase of Gaussian functions for modeling
%       x, y, and z coordinates of the cardiac dipole
%   delta_tetai: percentage of Gaussian center deviations added per beat
%   teta0: initial phase of the synthetic dipole
%
% Outputs:
%   dipole: structure containing the x, y, and z coordinates of the cardiac dipole
%   teta: vector containing the dipole phase
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

w = 2*pi*f; % Angular frequency
dt = 1/fs; % Time step

teta = zeros(1,N); % Initialize dipole phase vector
X = zeros(1,N); % Initialize x-coordinate vector
Y = zeros(1,N); % Initialize y-coordinate vector
Z = zeros(1,N); % Initialize z-coordinate vector

teta(1) = teta0; % Set initial dipole phase
d_alphai.x = alphai.x; % Initialize dipole amplitude deviations
d_alphai.y = alphai.y;
d_alphai.z = alphai.z;

d_tetai.x = tetai.x; % Initialize dipole phase deviations
d_tetai.y = tetai.y;
d_tetai.z = tetai.z;

d_bi.x = bi.x; % Initialize Gaussian width deviations
d_bi.y = bi.y;
d_bi.z = bi.z;

for i = 1 : N-1
    dtetaix = mod(teta(i) - d_tetai.x + pi , 2*pi) - pi; % Phase difference for x-coordinate
    dtetaiy = mod(teta(i) - d_tetai.y + pi , 2*pi) - pi; % Phase difference for y-coordinate
    dtetaiz = mod(teta(i) - d_tetai.z + pi , 2*pi) - pi; % Phase difference for z-coordinate

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
    teta(i+1) = teta(i) + w*dt;
    if teta(i+1) > pi % Beat transition
        teta(i+1) = teta(i+1) - 2*pi;
        
        % Apply stochastic deviations to parameters
        d_alphai.x = alphai.x * (1 + (rand - 0.5) * delta_alphai);
        d_alphai.y = alphai.y * (1 + (rand - 0.5) * delta_alphai);
        d_alphai.z = alphai.z * (1 + (rand - 0.5) * delta_alphai);

        d_tetai.x = tetai.x * (1 + (rand - 0.5) * delta_tetai);
        d_tetai.y = tetai.y * (1 + (rand - 0.5) * delta_tetai);
        d_tetai.z = tetai.z * (1 + (rand - 0.5) * delta_tetai);

        d_bi.x = bi.x * max(0, (1 + (rand - 0.5) * delta_bi));
        d_bi.y = bi.y * max(0, (1 + (rand - 0.5) * delta_bi));
        d_bi.z = bi.z * max(0, (1 + (rand - 0.5) * delta_bi));
        
        % Update angular frequency with stochastic deviation
        w = 2 * pi * f * max(0, (1 + (rand - 0.5) * f_deviations));
    end
end

% Store the dipole coordinates in the output structure
dipole.x = X;
dipole.y = Y;
dipole.z = Z;
