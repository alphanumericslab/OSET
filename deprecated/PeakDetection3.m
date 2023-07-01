function [peaks, mn, r] = PeakDetection3(x, fs, h, th, fmax)
% Deprecated: PeakDetection3 is deprecated. Use peak_detection_matched_filter instead.
    warning('Deprecated: PeakDetection3 is deprecated. Use peak_detection_matched_filter instead.');
    [peaks, mn, r] = peak_detection_matched_filter(x, fs, h, th, fmax);
end