function [DIP teta]= DipoleGeneratorAbnormal3(N,fs,rr,alphai,bi,tetai,teta0,STM,S0)
%
% [DIP teta]= DipoleGeneratorAbnormal3(N,fs,rr,alphai,bi,tetai,teta0,STM,S0)
% Synthetic cardiac dipole generator using the 'differential form' of the
% dipole equations. Refer to references of the toolbox for further details.
%
% inputs:
% N: signal length
% fs: sampling rate
% rr: rr interval time series
% alphai: structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% bi: structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% tetai: structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% teta0: initial phase of the synthetic dipole
% STM: the State Transition Matrix
% S0: initial state
% 
% Notes:
%     - For each entry of STM, Sij represents the probability of
%       going from state i to state j in the next beat
%     - Each row of STM should sum up to 1
%     - STM is usually asymmetric

%
% output:
% DIP: structure contaning the x, y, and z coordinates of the cardiac dipole
% teta: vector containing the dipole phase
%
%
% Open Source ECG Toolbox, version 2.0, April 2008
% Released under the GNU General Public License
% Copyright (C) 2008  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% Last update Gari Clifford April 16th 2006
%

if nargin < 9
    S0 = 1; % Start with a normal beat
end

if nargin < 8
 STM = 1;
 %STM = [0 1 ; 1 0]         % exact alternation: the T-wave will alternate in each beat
 STM1 = [.2 .8 ; .8 .2];   % probabilistic alternation with high probability of state transition
 STM2 = [.9 .1 ; .99 .01];
end

% make a time series of rr intervals resampled to the current sampling
% frequency
f=1./resample(rr,N,length(rr));
%f=mean(rr)*ones(1,N);

w = 2*pi*f;
dt = 1/fs;

%figure;plot(f);title('f')
%figure;plot(w);title('w')
%pause

teta = zeros(1,N);
X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);

CSTM = cumsum(STM,2);
L = size(STM,1);

teta(1) = teta0;
state = S0;
for i = 1:N-1;
    teta(i+1) = teta(i) + w(i)*dt;
    if(teta(i+1)>pi) % beat transition
        teta(i+1) = teta(i+1) - 2*pi;

    % Calculate heart rate-related probability transition matrix
    strt = i-(fs*10); if strt<1; strt=1; end
    hr = nanmean(f(strt):f(i));
    p=tanh((hr-95)/5);if p<0; p=0; end;
    STM = [1-p  p ; p 1-p];  
    CSTM = cumsum(STM,2);
    
        % new state estimation
        a = rand;
        if(a < CSTM(state,1))
            state = 1;
        else
            for j = 2:L,
                if(a >= CSTM(state,j-1) && a < CSTM(state,j))
                    state = j;
                    break;
                end
            end
        end
    end

    
    % simulate QT hysteresis
    strt = i-(fs*10); if strt<1; strt=1; end
    hr = nanmean(f(strt):f(i));
    qtc = sqrt(hr);
    qtc=1;
    
    dtetaix = mod(teta(i) - tetai(state).x + pi , 2*pi) - pi;
    dtetaiy = mod(teta(i) - tetai(state).y + pi , 2*pi) - pi;
    dtetaiz = mod(teta(i) - tetai(state).z + pi , 2*pi) - pi;

    if(i==1),
        X(i) = sum(alphai(state).x .* exp(-dtetaix .^2 ./ (2*bi(state).x .^ 2)));
        Y(i) = sum(alphai(state).y .* exp(-dtetaiy .^2 ./ (2*bi(state).y .^ 2)));
        Z(i) = sum(alphai(state).z .* exp(-dtetaiz .^2 ./ (2*bi(state).z .^ 2)));
    end

    X(i+1) = X(i) - dt*sum(w(i)*qtc*alphai(state).x ./ (bi(state).x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* bi(state).x .^ 2)));   % x state variable
    Y(i+1) = Y(i) - dt*sum(w(i)*qtc*alphai(state).y ./ (bi(state).y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* bi(state).y .^ 2)));   % y state variable
    Z(i+1) = Z(i) - dt*sum(w(i)*qtc*alphai(state).z ./ (bi(state).z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* bi(state).z .^ 2)));   % z state variable

end

DIP.x = X;
DIP.y = Y;
DIP.z = Z;