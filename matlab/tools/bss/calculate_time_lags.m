function [T0, T1] = calculate_time_lags(peaks, phase)
% calculate_time_lags - Calculate time lags for the pseudo-periodic component
% analysis algorithm (PiCA) using the R-peaks and the ECG phase signal.
%
% Usage:
%   [T0, T1] = calculate_time_lags(peaks, phase)
%
% Inputs:
%   peaks: Vector containing the indices of detected peaks in the input signal
%   phase: Vector representing the phase values associated with the peaks
%   (see phase_calculator)
%
% Outputs:
%   T0: Vector containing the time lags for the input peaks
%   T1: Vector containing the corresponding lagged time points
%
% Reference:
%   R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
%       electrocardiogram decomposition using periodic component analysis. IEEE
%       Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version CalculateTimeLags
% 
% Reza Sameni, 2008-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% PM time calculation
J = find(peaks);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1) - min(n1);

T1 = zeros(length(peaks) - prd - wlen, 1);
NN = length(T1);
for t = 1:NN
    df = abs(phase(t) - phase(t + prd - wlen : t + prd + wlen));
    [~, I] = min(df);
    T1(t) = t + prd + I - wlen - 1;
end

T0 = 1:NN;
