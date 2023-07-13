function [T0, T1] = synchronous_phase_samples(peaks)
% synchronous_phase_samples - Calculation of synchronous time instants from beat to
%   beat given a set of R-peaks, required for the periodic component analysis
%   (PiCA) algorithm.
%
% Usage:
%   [T0, T1] = synchronous_phase_samples(peaks)
%
% Inputs:
%   peaks: Vector of R-peak pulse train.
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
%   2023: Renamed from deprecated version SynchPhaseTimes2
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

I = find(peaks);
D = I(2:end) - I(1:end-1);

if length(I) < 3
    T0 = [];
    T1 = [];
else
    start = I(1);
    stop = I(end-1);

    T1 = zeros(1, stop - start + 1);
    k = 1;
    for t = start:stop
        T1(t - start + 1) = I(k + 1) + round((t - I(k)) * D(k + 1) / D(k));
        if (t >= I(k + 1))
            k = k + 1;
        end
    end
    T0 = start:stop;
end
