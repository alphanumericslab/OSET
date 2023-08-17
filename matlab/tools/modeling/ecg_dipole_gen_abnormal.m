function [dipole, teta]= ecg_dipole_gen_abnormal(N,fs,f,alphai,bi,tetai,teta0,STM,S0)
% ecg_dipole_gen_abnormal - Synthetic cardiac dipole generator using the
%   state-space form of the sum of Gaussian dipole model. Refer to references
%   for further details.
%
% Usage:
%   [dipole teta]= ecg_dipole_gen_abnormal(N,fs,f,alphai,bi,tetai,teta0,STM,S0)
% 
% inputs:
%   N: signal length
%   fs: sampling rate
%   f: average heart rate (Hz)
%   alphai: vector structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   bi: vector structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   tetai: vector structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   teta0: vector initial phase of the synthetic dipole
%   STM: the state transition matrix from one beat type to another
%   S0: initial state
% 
% output:
%   dipole: structure contaning the x, y, and z coordinates of the cardiac dipole
%   teta: vector containing the dipole phase
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

teta = zeros(1,N); % Initialize dipole phase vector
X = zeros(1,N); % Initialize x-coordinate vector
Y = zeros(1,N); % Initialize y-coordinate vector
Z = zeros(1,N); % Initialize z-coordinate vector

CSTM = cumsum(STM,2); % Cumulative sum of state transition matrix
L = size(STM,1); % Number of states

teta(1) = teta0; % Set initial dipole phase
state = S0; % Set initial state
for i = 1:N-1
    teta(i+1) = teta(i) + w*dt; % Update dipole phase
    
    % Check for beat transition
    if teta(i+1) > pi
        teta(i+1) = teta(i+1) - 2*pi;
        
        % Estimate new state using random transition
        a = rand;
        if a < CSTM(state,1)
            state = 1;
        else
            for j = 2:L
                if a >= CSTM(state,j-1) && a < CSTM(state,j)
                    state = j;
                    break;
                end
            end
        end
    end
    
    dtetaix = mod(teta(i) - tetai(state).x + pi , 2*pi) - pi; % Phase difference for x-coordinate
    dtetaiy = mod(teta(i) - tetai(state).y + pi , 2*pi) - pi; % Phase difference for y-coordinate
    dtetaiz = mod(teta(i) - tetai(state).z + pi , 2*pi) - pi; % Phase difference for z-coordinate

    if i == 1
        % Compute initial values for the first iteration
        X(i) = sum(alphai(state).x .* exp(-dtetaix .^2 ./ (2*bi(state).x .^ 2)));
        Y(i) = sum(alphai(state).y .* exp(-dtetaiy .^2 ./ (2*bi(state).y .^ 2)));
        Z(i) = sum(alphai(state).z .* exp(-dtetaiz .^2 ./ (2*bi(state).z .^ 2)));
    end

    % Update x, y, and z state variables
    X(i+1) = X(i) - dt*sum(w*alphai(state).x ./ (bi(state).x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* bi(state).x .^ 2)));
    Y(i+1) = Y(i) - dt*sum(w*alphai(state).y ./ (bi(state).y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* bi(state).y .^ 2)));
    Z(i+1) = Z(i) - dt*sum(w*alphai(state).z ./ (bi(state).z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* bi(state).z .^ 2)));

end

% Store the dipole coordinates in the output structure
dipole.x = X;
dipole.y = Y;
dipole.z = Z;
