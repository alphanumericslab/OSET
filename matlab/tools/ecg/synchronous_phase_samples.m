function [T0, T1] = synchronous_phase_samples(peaks, varargin)
% synchronous_phase_samples - Calculation of synchronous time instants from beat to
%   beat given a set of R-peaks, required for the periodic component analysis
%   (PiCA) algorithm.
%
% Usage:
%   [T0, T1] = synchronous_phase_samples(peaks, phase)
%
% Inputs:
%   peaks: Vector of R-peak pulse train.
%   phase (optional): The ECG phase obtained from phase_calculator
%
% Outputs:
%   T0: First (or reference) time instant vector.
%   T1: Second time vector having synchronous phases with T0.
%
% Reference and usage:
% - R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel electrocardiogram
%   decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering,
%   55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2009: First release
%   2023: Combined deprecated versions SynchPhaseTimes and SynchPhaseTimes2
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

I_peaks = find(peaks);
D = diff(I_peaks);
if nargin > 1 && ~isempty(varargin{1}) % Use pre-calculated ECG phase. Identical to deprecated function SynchPhaseTimes
    phase = varargin{1};

    L = length(peaks);
    prd = round(mean(D));
    wlen = max(D)-min(D);

    T1 = zeros(1,L-prd-wlen);
    NN = length(T1);
    for t = 1:NN
        df = abs(phase(t) - phase(max(t+prd-wlen,1):min(t+prd+wlen,L)));
        [~, I] = min(df);
        T1(t) = t + prd + I - wlen -1;
    end
    T1 = max(T1,1);
    T1 = min(T1,NN);
    T0 = 1 : NN;

else % identical to deprecated function SynchPhaseTimes2
    if length(I_peaks) < 3
        T0 = [];
        T1 = [];
    else
        start = I_peaks(1);
        stop = I_peaks(end-1);

        T1 = zeros(1, stop - start + 1);
        k = 1;
        for t = start : stop
            T1(t - start + 1) = I_peaks(k + 1) + round((t - I_peaks(k)) * D(k + 1) / D(k));
            if t >= I_peaks(k + 1)
                k = k + 1;
            end
        end
        T0 = start : stop;
    end
end