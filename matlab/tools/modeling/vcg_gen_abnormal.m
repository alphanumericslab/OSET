function [vcg, phi]= vcg_gen_abnormal(N,fs,f,alpha,b,theta,theta0,STM,S0)
%
% vcg_gen_abnormal - Synthetic cardiac dipole/vectorcardiogram(VCG) generator using the
%   state-space form of the sum of Gaussian dipole/VCG model. Refer to references
%   for further details.
%
% Usage:
%   [vcg, phi]= vcg_gen_abnormal(N,fs,f,alpha,b,theta,theta0,STM,S0)
%
% inputs:
%   N: signal length
%   fs: sampling rate
%   f: average heart rate (Hz)
%   alpha: vector structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole/VCG
%   b: vector structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole/VCG
%   theta: vector structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole/VCG
%   theta0: vector initial phase of the synthetic dipole/VCG
%   STM: the state transition matrix from one beat type to another
%   S0: initial state
%
% output:
%   vcg: structure contaning the x, y, and z coordinates of the cardiac dipole/VCG
%   phi: vector containing the dipole/VCG phase
%
% Notes:
%     - For each entry of STM, Sij represents the probability of
%       going from state i to state j in the next beat
%     - Each row of STM should sum up to 1
%     - STM is usually asymmetric
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
%   2008: First release
%   2023: Renamed from the deprecated version DipoleGeneratorAbnormal

w = 2*pi*f; % Angular frequency
dt = 1/fs; % Time step

phi = zeros(1,N); % Initialize dipole/VCG phase vector
X = zeros(1,N); % Initialize x-coordinate vector
Y = zeros(1,N); % Initialize y-coordinate vector
Z = zeros(1,N); % Initialize z-coordinate vector

% CSTM = cumsum(STM,2) % Cumulative sum of state transition matrix
% L = size(STM,1); % Number of states

phi(1) = theta0; % Set initial dipole/VCG phase
state = S0; % Set initial state
for i = 1:N-1
    phi(i+1) = phi(i) + w*dt; % Update dipole/VCG phase

    % Check for beat transition
    if phi(i+1) > pi
        phi(i+1) = phi(i+1) - 2*pi;

        % Update the state according to the provided STM
        state = markov_model_next_state_gen(STM, state);

    end

    dtetaix = mod(phi(i) - theta(state).x + pi , 2*pi) - pi; % Phase difference for x-coordinate
    dtetaiy = mod(phi(i) - theta(state).y + pi , 2*pi) - pi; % Phase difference for y-coordinate
    dtetaiz = mod(phi(i) - theta(state).z + pi , 2*pi) - pi; % Phase difference for z-coordinate

    if i == 1
        % Compute initial values for the first iteration
        X(i) = sum(alpha(state).x .* exp(-dtetaix .^2 ./ (2*b(state).x .^ 2)));
        Y(i) = sum(alpha(state).y .* exp(-dtetaiy .^2 ./ (2*b(state).y .^ 2)));
        Z(i) = sum(alpha(state).z .* exp(-dtetaiz .^2 ./ (2*b(state).z .^ 2)));
    end

    % Update x, y, and z state variables
    X(i+1) = X(i) - dt*sum(w*alpha(state).x ./ (b(state).x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* b(state).x .^ 2)));
    Y(i+1) = Y(i) - dt*sum(w*alpha(state).y ./ (b(state).y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* b(state).y .^ 2)));
    Z(i+1) = Z(i) - dt*sum(w*alpha(state).z ./ (b(state).z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* b(state).z .^ 2)));

end

% Store the dipole/VCG coordinates in the output structure
vcg.x = X;
vcg.y = Y;
vcg.z = Z;

end
