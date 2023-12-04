function [peaks, mn, r] = peak_det_matched_filter(x, fs, h, th, fmax)
% peak_det_matched_filter - R-peak detector based on a matched filter
%
% Syntax: [peaks, mn, r] = peak_det_matched_filter(x, fs, h, th, fmax)
%
% Inputs:
%       x:      Vector of input data
%       fs:     Sampling rate
%       h:      Template waveform
%       th:     Detection threshold
%       fmax:   Maximum expected frequency of the R-peaks
%
% Outputs:
%       peaks:  Vector of R-peak impulse train
%       mn:     Mean waveform of detected R-peaks
%       r:      Filtered output after matched filtering
%
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version PeakDetection3()
%
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

N = length(x);
L = length(h);

h = h(end:-1:1); % Reverse the template waveform

w = round(L/2);

r = filter(h, 1, [x zeros(1, w-1)]); % Matched filtering
r = r(w : N+w-1); % Trim the filtered output to match the length of x
r(r < th*max(r)) = 0; % Thresholding: set values below th*max(r) to 0

peaks = zeros(1, length(x));

if nargout == 1 % only peaks output requested
    wlen2 = round(fs/fmax);
    I = find(r > 0);

    for i = 1 : length(I)
        ind = max(1, I(i)-wlen2) : min(N, I(i)+wlen2);
        if max(r(ind)) == r(I(i))
            peaks(I(i)) = 1;
        end
    end
else % peaks and mn outputs requested
    wlen2 = round(fs/fmax);
    I = find(r > 0);

    seg = zeros(length(I), L);
    for i = 1 : length(I)
        ind = max(1, I(i)-wlen2) : min(N, I(i)+wlen2);
        if max(r(ind)) == r(I(i))
            peaks(I(i)) = 1;
            sg = x(max(1, I(i)-w+1) : min(N, I(i)+w));
            seg(i, :) = [sg zeros(1, L-length(sg))];
        end
    end

    mn0 = mean(seg, 1);

    % Robust weighted averaging
    noise = seg - mn0(ones(size(seg, 1), 1), :);
    vr = var(noise, [], 2);
    sm = sum(1 ./ vr);
    weight = 1 ./ (vr * sm);
    mn = weight' * seg;

    mn = sqrt(sum(h.^2)) * mn / sqrt(sum(mn.^2));
end
