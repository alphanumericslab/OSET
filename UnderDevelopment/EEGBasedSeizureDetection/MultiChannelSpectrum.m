function H = MultiChannelSpectrum(x, fs, nfft, h)

H = dspdata.psd;
H = H(ones(1, size(x,1)));

for i = 1:size(x, 1),
    H(i) = psd(h, x(i,:), 'Fs', fs, 'NFFT', nfft - 1, 'SpectrumType','onesided'); % Calculate the PSD 
end
