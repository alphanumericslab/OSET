function [DIP, teta]= DipoleGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
%
% [DIP teta]= DipoleGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
% Synthetic cardiac dipole generator with beat-wise stochatic deviations.
%
% inputs:
% N: signal length
% fs: sampling rate
% f: average heart rate (Hz)
% f_deviations: percentage of beat-wise heart rate deviations (Hz)
% alphai: structure contaning the amplitudes of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% delta_alphai: the percentage of amplitude deviations added per beat
% bi: structure contaning the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% delta_bi: the percentage of gaussian wave width deviations added per beat
% tetai: structure contaning the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
% delta_tetai: the percentage of gaussian center deviations added per beat
% teta0: initial phase of the synthetic dipole
%

%
% output:
% DIP: structure contaning the x, y, and z coordinates of the cardiac dipole
% teta: vector containing the dipole phase
%
%
% The Open Source Electrophysiological Toolbox (OSET), version 3.14, April 2022
% URL: https://github.com/alphanumericslab/OSET
% Copyright (C) 2022  Reza Sameni
% reza.sameni@gmail.com

w = 2*pi*f;
dt = 1/fs;

teta = zeros(1,N);
X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);

teta(1) = teta0;
d_alphai.x = alphai.x;
d_alphai.y = alphai.y;
d_alphai.z = alphai.z;

d_tetai.x = tetai.x;
d_tetai.y = tetai.y;
d_tetai.z = tetai.z;

d_bi.x = bi.x;
d_bi.y = bi.y;
d_bi.z = bi.z;
for i = 1 : N-1
    dtetaix = mod(teta(i) - d_tetai.x + pi , 2*pi) - pi;
    dtetaiy = mod(teta(i) - d_tetai.y + pi , 2*pi) - pi;
    dtetaiz = mod(teta(i) - d_tetai.z + pi , 2*pi) - pi;

    if(i==1)
        X(i) = sum(d_alphai.x .* exp(-dtetaix .^2 ./ (2*d_bi.x .^ 2)));
        Y(i) = sum(d_alphai.y .* exp(-dtetaiy .^2 ./ (2*d_bi.y .^ 2)));
        Z(i) = sum(d_alphai.z .* exp(-dtetaiz .^2 ./ (2*d_bi.z .^ 2)));
    end

    X(i+1) = X(i) - dt*sum(w * d_alphai.x ./ (d_bi.x .^ 2) .* dtetaix .* exp(-dtetaix .^2 ./ (2* d_bi.x .^ 2)));   % x state variable
    Y(i+1) = Y(i) - dt*sum(w * d_alphai.y ./ (d_bi.y .^ 2) .* dtetaiy .* exp(-dtetaiy .^2 ./ (2* d_bi.y .^ 2)));   % y state variable
    Z(i+1) = Z(i) - dt*sum(w * d_alphai.z ./ (d_bi.z .^ 2) .* dtetaiz .* exp(-dtetaiz .^2 ./ (2* d_bi.z .^ 2)));   % z state variable

    % Next beat
    teta(i+1) = teta(i) + w*dt;
    if(teta(i+1) > pi) % beat transitions
        teta(i+1) = teta(i+1) - 2*pi;
        d_alphai.x = alphai.x * (1 + (rand - 0.5) * delta_alphai);
        d_alphai.y = alphai.y * (1 + (rand - 0.5) * delta_alphai);
        d_alphai.z = alphai.z * (1 + (rand - 0.5) * delta_alphai);

        d_tetai.x = tetai.x * (1 + (rand - 0.5) * delta_tetai);
        d_tetai.y = tetai.y * (1 + (rand - 0.5) * delta_tetai);
        d_tetai.z = tetai.z * (1 + (rand - 0.5) * delta_tetai);

        d_bi.x = bi.x * max(0, (1 + (rand - 0.5) * delta_bi));
        d_bi.y = bi.y * max(0, (1 + (rand - 0.5) * delta_bi));
        d_bi.z = bi.z * max(0, (1 + (rand - 0.5) * delta_bi));
        w = 2 * pi * f * max(0, (1 + (rand - 0.5) * f_deviations));
    end
end

DIP.x = X;
DIP.y = Y;
DIP.z = Z;