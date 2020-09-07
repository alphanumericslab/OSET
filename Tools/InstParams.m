function [A, f0, bw, fskew, skew, kurt, f1] = InstParams(signal, fs, wlen, window, th, flag)
% 
% This function calculates instantaneous time/frequency domain features of an
% input signal in causal and non-causal modes
% 
% [A, f0, bw, fskew, skew, kurt, f1] = InstParams(signal, fs, wlen, window, th, flag)
% Inputs:
%   signal: input signal as row vector
%   fs: sampling frequency
%   wlen: sliding window length used for parameter calculation
%   window: a row-wise vector of size wlen used for signal windowing before
%           parameter calculation
%   th: instantaneous amplitude threshold for frequency domain parameter
%           calculation to avoid noisy reports when analytical signal envelope 
%           is too weak. Set to 0 to calculate in all amplitudes.
%   flag: Set to 0 for causal and 1 for non-causal calculation
%
% Outputs:
%   A: Analytical form signal envelope
%   f0: instantaneous frequency using energy density weighting
%   bw: instantaneous bandwidth using energy density weighting
%   fskew: instantaneous frequency skewness using energy density weighting
%   skew: instantaneous time-domain skewness using energy density weighting
%   kurt: instantaneous time-domain kurtosis using energy density weighting
%   f1: instantaneous frequency using Hilbert transform

%
% Reza Sameni, March 2015
% Modified:
%   Sep 2020: Help added
% 
% reza.sameni@gmail.com
% Open-Source Electrophysiological Toolbox (OSET)
% www.oset.ir
% Git repository: https://gitlab.com/rsameni/OSET/

if(isempty(flag))
    flag = 0; % causal
end

% make the signals row-wise
signal = signal(:)';
window = window(:)';

A = zeros(size(signal));
f0 = zeros(size(signal));

h = hilbert(signal);
f1 = fs*diff(unwrap(angle(h)))/(2*pi);
f1 = [f1(1) f1];

bw = zeros(size(signal));
fskew = zeros(size(signal));
skew = zeros(size(signal));
kurt = zeros(size(signal)) + 3; % 3 for gaussian
N = length(signal);        % signal length
if(flag == 0) % causal
    for i = wlen : N
        x = signal(i-wlen+1 : i);
        x = window .* x; % windowing
        NFFT = length(x);
        halfNFFT = round(NFFT/2);
        X = fft(x, NFFT);
        XX = abs(X(1:halfNFFT)).^2;
        X_energy_density = XX / mean(XX);
        f = (fs/2)*(0:halfNFFT-1)/halfNFFT;
        A(i) = sqrt(mean(abs(x).^2));
        skew(i) = skewness(x);
        kurt(i) = kurtosis(x);
        if A(i) > th
            f0(i) = mean(f .* X_energy_density);
            bw(i) = sqrt(mean( (f - f0(i)).^2 .* X_energy_density));
            fskew(i) = nthroot(mean( (f - f0(i)).^3 .* X_energy_density ), 3);
        end
    end
    
elseif(flag == 1)% non-causal
    for i = ceil(wlen/2) : N-floor(wlen/2)
        x = signal(i-ceil(wlen/2)+1 : i+floor(wlen/2));
        x = window .* x; % windowing
        NFFT = length(x);
        halfNFFT = round(NFFT/2);
        X = fft(x,NFFT);
        XX = abs(X(1:halfNFFT)).^2;
        X_energy_density = XX / mean(XX);
        f = (fs/2)*(0:halfNFFT-1)/halfNFFT;
        A(i) = sqrt(mean(abs(x).^2));
        skew(i) = skewness(x);
        kurt(i) = kurtosis(x);
        if A(i) > th
            f0(i) = mean(f .* X_energy_density);
            bw(i) = sqrt(mean( (f - f0(i)).^2 .* X_energy_density ));
            fskew(i) = power(mean( (f - f0(i)).^3 .* X_energy_density ), 1/3);
        end
    end
    
end
