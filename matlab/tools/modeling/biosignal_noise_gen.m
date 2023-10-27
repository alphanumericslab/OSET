function noise = biosignal_noise_gen(noisetype, SignalPower, snr, N, varargin)
% biosignal_noise_gen - Generate biosignal noise based on various types.
%
% Description:
%   This function generates biosignal noise based on different modes, such
%   as white noise, colored noise, real muscle artifacts, real electrode
%   movements, real baseline wander, or a mixture of real baseline wander,
%   electrode movements, muscle artifacts and white noise. The real noise
%   samples have been adopted from The MIT-BIH Noise Stress Test Database: 
%   https://physionet.org/content/nstdb/1.0.0/
% 
% Usage:
%       WN =  biosignal_noise_gen(0,SignalPower,snr,N,seed);
%       CN =  biosignal_noise_gen(1,SignalPower,snr,N,fs,beta,seed);
%       MA =  biosignal_noise_gen(2,SignalPower,snr,N,fs,seed);
%       EM =  biosignal_noise_gen(3,SignalPower,snr,N,fs,seed);
%       BW =  biosignal_noise_gen(4,SignalPower,snr,N,fs,seed);
%       MX =  biosignal_noise_gen(5,SignalPower,snr,N,fs,[w_bw,w_em,w_ma],seed);
%
% Modes:
%   0: White Noise (WN)
%   1: Colored Noise (CN)
%   2: Real Muscle Artifacts (MA)
%   3: Real Electrode Movements (EM)
%   4: Real Baseline Wander (BW)
%   5: Mixture of Real BW, EM, and MA (MX)
%
% Inputs:
%   noisetype: Type of noise to generate (0-5).
%   SignalPower: Desired signal power (set to mean(x.^2) for data vector x).
%   snr: Desired Signal-to-Noise Ratio (DNR) in dB.
%   N: Number of samples.
%   Additional input arguments based on the mode.
%       For mode 0: seed (optional). random seed
%       For mode 1: fs (sampling frequency), beta (noise coloring factor), seed (optional)
%       For mode 2, 3, 4: fs (sampling frequency), seed (optional)
%       For mode 5: fs (sampling frequency), [w_bw, w_em, w_ma] (weights of BW, EM, MA), seed (optional)
%
% Output:
%   noise: Column vector of generated noise.
%
% References:
%   Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%   Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%   ECG Signals. EURASIP Journal on Advances in Signal Processing, 2007(1).
% 
% Revision History:
%   2006: First release.
%   2023: Documented and renamed from deprecated version NoiseGenerator.
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

switch noisetype
    case 0     % White Noise
        if nargin == 5
            rng(varargin{1});
        end
        NoisePower = SignalPower / 10^(snr/10);
        noise = sqrt(NoisePower) * randn(N, 1);

    case 1     % Colored Noise
        fs = varargin{1};
        beta = varargin{2};
        if nargin == 7
            rng(varargin{3});
        end
        NoisePower = SignalPower / 10^(snr/10);
        noise = ColoredNoise(sqrt(NoisePower), N, fs, beta);

    case 2     % Real Muscle Artifacts
        fs = varargin{1};
        NoisePower = SignalPower / 10^(snr/10);
        load('MA.mat', 'MA');
        artifact = MA(:, 2);
        if fs ~= 360.0
            artifact = resample(artifact, fs, 360);
        end
        M = length(artifact);
        if nargin == 6
            n0 = max(1, min(varargin{2}, length(artifact)-N+1));
        else
            n0 = 1;
        end
        if M-n0 < N
            artifact = wextend('1','sym',artifact,N-(M-n0),'r');
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower) * (artifact - mean(artifact)) / std(artifact, 1);

    case 3     % Real Electrode Movements
        fs = varargin{1};
        NoisePower = SignalPower / 10^(snr/10);
        load('EM.mat', 'EM');
        artifact = EM(:, 3);
        if fs ~= 360.0
            artifact = resample(artifact, fs, 360);
        end
        M = length(artifact);
        if nargin == 6
            n0 = max(1, min(varargin{2}, length(artifact)-N+1));
        else
            n0 = 1;
        end
        if M-n0 < N
            artifact = wextend('1','sym',artifact,N-(M-n0),'r');
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower) * (artifact - mean(artifact)) / std(artifact, 1);

    case 4     % Real Baseline Wander
        fs = varargin{1};
        NoisePower = SignalPower / 10^(snr/10);
        load('BW.mat', 'BW');
        artifact = BW(:, 3);
        if fs ~= 360.0
            artifact = resample(artifact, fs, 360);
        end
        M = length(artifact);
        if nargin == 6
            n0 = max(1, min(varargin{2}, length(artifact)-N+1));
        else
            n0 = 1;
        end
        if M-n0 < N
            artifact = wextend('1','sym',artifact,N-(M-n0),'r');
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower) * (artifact - mean(artifact)) / std(artifact, 1);

    case 5     % Mixture of Real BW, EM, and MA
        fs = varargin{1};
        w = varargin{2};
        w_bw = w(1);       % Weight of baseline wander noise in the generated noise
        w_em = w(2);       % Weight of electrode movement noise in the generated noise
        w_ma = w(3);       % Weight of muscle artifact noise in the generated noise
        NoisePower = SignalPower / 10^(snr/10);
        load('BW.mat', 'BW'); bw = BW(:, 3);    bw = (bw - mean(bw)) / std(bw);
        load('EM.mat', 'EM'); em = EM(:, 3);    em = (em - mean(em)) / std(em);
        load('MA.mat', 'MA'); ma = MA(:, 3);    ma = (ma - mean(ma)) / std(ma);
        artifact = (w_bw * bw + w_em * em + w_ma * ma) / (w_bw + w_em + w_ma);
        artifact = resample(artifact, fs, 360);
        M = length(artifact);
        if nargin == 7
            n0 = max(1, min(varargin{3}, length(artifact)-N+1));
        else
            n0 = 1;
        end
        if M-n0 < N
            artifact = wextend('1','sym',artifact,N-(M-n0),'r');
        end
        artifact = artifact(n0:N+n0-1)';
        noise = sqrt(NoisePower) * (artifact - mean(artifact)) / std(artifact, 1);
end

noise = noise(:);
