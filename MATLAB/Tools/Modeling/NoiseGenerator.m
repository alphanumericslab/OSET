function noise =  NoiseGenerator(noisetype,SignalPower,SNR,N,varargin);
%
% noise =  NoiseGenerator(noisetype,SignalPower,SNR,N,Param1,Param2,...),
% ECG noise generator
%
% Usage:
%       WN =  NoiseGenerator(0,SignalPower,SNR,N,seed);
%       CN =  NoiseGenerator(1,SignalPower,SNR,N,fs,beta,seed);
%       MA =  NoiseGenerator(2,SignalPower,SNR,N,fs,seed);
%       EM =  NoiseGenerator(3,SignalPower,SNR,N,fs,seed);
%       BW =  NoiseGenerator(4,SignalPower,SNR,N,fs,seed);
%       MX =  NoiseGenerator(5,SignalPower,SNR,N,fs,[w_bw,w_em,w_ma],seed);
%
% inputs:
% noisetype
%       0:     white noise (WN)
%       1:     colored noise (CN)
%       2:     real muscle artifacts (MA)
%       3:     real electrode movements (EM)
%       4:     real baseline wander (BW)
%       5:     mixture of real baseline wander, electrode movements, muscle artifacts (MX)
% SignalPower: The desired signal power. set to mean(x.^2) for the data vector x
% SNR: The desired SNR
% N: Number of samples
% fs: Sampling frequency required for noisetype = 1,...,5
% beta: Noise coloring factor required for noisetype = 1. beta = 0 (white noise),
%       beta = 1 (pink noise), beta = 2 (brown noise or random walk)
% seed(optional): Random seed for the noise vector. For noisetype =
%       2,...,5 seed is the initial random starting point in the real recorded
%       noises
% [w_bw,w_em,w_ma](optional): The weighting factors of BW, EM, and MA noise (only for noisetype = 5).
%
% output:
% noise: Column vector of noise
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.

switch noisetype
    case 0     % white noise
        if (nargin==5),
            randn('seed',varargin{1});
        end
        NoisePower = SignalPower / 10^(SNR/10);
        noise = sqrt(NoisePower)*randn(N,1);

    case 1     % colored noise
        fs = varargin{1};
        beta = varargin{2};
        if (nargin==7),
            randn('seed',varargin{3});
        end
        NoisePower = SignalPower / 10^(SNR/10);
        noise = ColoredNoise(sqrt(NoisePower),N,fs,beta);

    case 2     % real muscle artifacts
        fs = varargin{1};
        NoisePower = SignalPower / 10^(SNR/10);
        load('MA.mat');artifact = MA(:,2);
        artifact = resample(artifact,fs,360);
        if (nargin==6),
            n0 = max(1,min(varargin{2},length(artifact)-N+1));
        else
            n0 = 1;
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);

    case 3     % real electrode movements
        fs = varargin{1};
        NoisePower = SignalPower / 10^(SNR/10);
        load('EM.mat');artifact = EM(:,3);
        artifact = resample(artifact,fs,360);
        if (nargin==6),
            n0 = max(1,min(varargin{2},length(artifact)-N+1));
        else
            n0 = 1;
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);

    case 4     % real baseline wander
        fs = varargin{1};
        NoisePower = SignalPower / 10^(SNR/10);
        load('BW.mat');artifact = BW(:,3);
        artifact = resample(artifact,fs,360);
        if (nargin==6),
            n0 = max(1,min(varargin{2},length(artifact)-N+1));
        else
            n0 = 1;
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);

    case 5     % mixture of real baseline wander, electrode movements, muscle artifacts
        fs = varargin{1};
        w = varargin{2};
        w_bw = w(1);       % weight of baseline wander noise in the generated noise
        w_em = w(2);       % weight of electrode movement noise in the generated noise
        w_ma = w(3);       % weight of muscle artifact noise in the generated noise
        NoisePower = SignalPower / 10^(SNR/10);
        load('BW.mat'); bw = BW(:,3);    bw = (bw-mean(bw))/std(bw);
        load('EM.mat'); em = EM(:,3);    em = (em-mean(em))/std(em);
        load('MA.mat'); ma = MA(:,3);    ma = (ma-mean(ma))/std(ma);
        artifact = (w_bw*bw + w_em*em + w_ma*ma)/(w_bw + w_em + w_ma);
        artifact = resample(artifact,fs,360);
        if (nargin==7),
            n0 = max(1,min(varargin{3},length(artifact)-N+1));
        else
            n0 = 1;
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower)*(artifact - mean(artifact))/std(artifact,1);
end

noise = noise(:);