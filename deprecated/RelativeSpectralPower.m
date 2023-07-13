function ratios = RelativeSpectralPower(x, fs, freqs, Q_factor)
% RelativeSpectralPower has been deprecated. Use relative_spectral_power instead.
warning('RelativeSpectralPower has been deprecated. Use relative_spectral_power instead.');
ratios = relative_spectral_power(x, fs, freqs, Q_factor);