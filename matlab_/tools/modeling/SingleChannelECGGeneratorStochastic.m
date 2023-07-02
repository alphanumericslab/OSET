function [ECG, teta]= SingleChannelECGGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
%
% [DIP teta]= SingleChannelECGGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
% Synthetic single-channel ECG generator with beat-wise stochatic deviations.
%
% inputs:
% N: signal length
% fs: sampling rate
% f: average heart rate (Hz)
% f_deviations: percentage of beat-wise heart rate deviations (Hz)
% alphai: amplitudes of Gaussian functions used for ECG modeling
% delta_alphai: the percentage of amplitude deviations added per beat
% bi: widths of Gaussian functions used for ECG modeling
% delta_bi: the percentage of gaussian wave width deviations added per beat
% tetai: phase of Gaussian functions used for ECG modelin
% delta_tetai: the percentage of gaussian center deviations added per beat
% teta0: initial phase of the synthetic dipole
%

%
% output:
% ECG: Single-channel ECG
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
ECG = zeros(1,N);

teta(1) = teta0;
d_alphai = alphai;
d_tetai = tetai;
d_bi = bi;
for i = 1 : N-1
    dtetai = mod(teta(i) - d_tetai + pi , 2*pi) - pi;

    if(i==1)
        ECG(i) = sum(d_alphai .* exp(-dtetai .^2 ./ (2*d_bi .^ 2)));
    end

    ECG(i+1) = ECG(i) - dt*sum(w * d_alphai ./ (d_bi .^ 2) .* dtetai .* exp(-dtetai .^2 ./ (2* d_bi .^ 2)));
    
    % Next beat
    teta(i+1) = teta(i) + w * dt;
    if(teta(i+1) > pi) % beat transitions
        teta(i+1) = teta(i+1) - 2*pi;
        d_alphai = alphai * (1 + (rand - 0.5) * delta_alphai);

        d_tetai = tetai * (1 + (rand - 0.5) * delta_tetai);

        d_bi = bi * max(0, (1 + (rand - 0.5) * delta_bi));
        w = 2 * pi * f * max(0, (1 + (rand - 0.5) * f_deviations));
    end
end