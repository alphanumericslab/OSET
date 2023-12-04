function noise = NoiseGeneratorSimple(noisetype, NoisePower, fs, N)
% 
% NoiseGeneratorSimple - Generate simplified noise for simulation purposes.
% Note:
%   NoiseGeneratorSimple has been deprecated. Use biosignal_noise_gen instead.'
% Syntax:
%   noise = NoiseGeneratorSimple(noisetype, NoisePower, fs, N)
%
% Inputs:
%   noisetype: Type of noise to generate ('WHITE', 'MA', 'EM')
%   NoisePower: Power of the generated noise
%   fs: Sampling frequency
%   N: Number of samples
%
% Outputs:
%   noise: Generated noise signal.
%
% Description:
%   This function generates simplified noise signals for simulation purposes.
%   It includes white noise, muscle artifacts (MA), and electrode movements (EM)
%
% References:
%   Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%   Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%   ECG Signals. EURASIP Journal on Advances in Signal Processing, 2007(1)
%
% Revision History:
%   2006: First release.
%   2023: Documented and renamed from deprecated version NoiseGenerator
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

warning('NoiseGeneratorSimple has been deprecated. Use biosignal_noise_gen instead.');

switch noisetype
    case 'WHITE'
        % Generate white noise
        noise = sqrt(NoisePower) * randn(1, N);
        
    case 'MA'
        % Generate muscle artifact noise
        load('MA.mat', 'MA');
        ma = MA(:, 3)';
        if fs ~= 360.0
            ma = resample(ma, fs, 360.0);
        end
        ma = ma - LPFilter(ma, 1.0 / fs);
        ma = ma(1:N);
        ma = (ma - mean(ma)) / std(ma);
        noise = sqrt(NoisePower) * ma;
        
    case 'EM'
        % Generate electrode movement noise
        load('EM.mat', 'EM');
        em = EM(:, 3)';
        if fs ~= 360.0
            em = resample(em, fs, 360.0);
        end
        em = em - LPFilter(em, 1.0 / fs);
        em = em(1:N);
        em = (em - mean(em)) / std(em);
        noise = sqrt(NoisePower) * em;
end
