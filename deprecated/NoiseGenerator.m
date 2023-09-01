function noise = NoiseGenerator(noisetype, SignalPower, SNR, N, varargin)
% NoiseGenerator has been deprecated. Use biosignal_noise_gen instead.
warning('NoiseGenerator has been deprecated. Use biosignal_noise_gen instead.');
noise = biosignal_noise_gen(noisetype, SignalPower, SNR, N, varargin{:});