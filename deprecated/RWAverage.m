function [mn, vr_mn, md, vr_md] = RWAverage(x)
% RWAverage has been deprecated. Use robust_weighted_average instead.
warning('RWAverage has been deprecated. Use robust_weighted_average instead.');
[mn, vr_mn, md, vr_md] = robust_weighted_average(x);