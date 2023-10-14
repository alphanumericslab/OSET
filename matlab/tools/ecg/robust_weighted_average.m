function [mn, vr_mn, md, vr_md] = robust_weighted_average(x)
% robust_weighted_average - Robust weighted averaging
%
% Usage:
%   [mn, vr_mn, md, vr_md] = robust_weighted_average(x)
%
% Inputs:
%   x: An (N x T) matrix containing N ensembles of a noisy event-related signal of length T
%
% Outputs:
%   mn: The robust weighted average over the N rows of x
%   vr_mn: The variance of the average beat across the N rows of x
%   md: The robust weighted median over the N rows of x
%   vr_md: The variance of the median beat across the N rows of x
%
% Reference:
%   J.M. Leski. Robust weighted averaging of biomedical signals. IEEE
%       Trans. Biomed. Eng., 49(8):796-804, 2002.
%
% Revision History:
%   2008: First release
%   2019: Return the average beat variance
%   2021: Return the median beat and its variance
%   2023: Renamed from deprecated version RWAverage
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

if (size(x, 1) > 1)
    % Average beat (mean)
    mn0 = mean(x, 1);
    noise0 = x - mn0(ones(size(x, 1), 1), :);
    vr = var(noise0, [], 2);
    sm = sum(1 ./ vr);
    weight = 1 ./ (vr * sm);
    mn = weight' * x;
    noise = x - mn(ones(size(x, 1), 1), :);
    vr_mn = var(noise, [], 1);

    % Average beat (median)
    md0 = median(x, 1);
    noise0 = x - md0(ones(size(x, 1), 1), :);
    vr = var(noise0, [], 2);
    sm = sum(1 ./ vr);
    weight = 1 ./ (vr * sm);
    md = weight' * x;
    noise = x - md(ones(size(x, 1), 1), :);
    vr_md = var(noise, [], 1);
else
    mn = x;
    md = x;
    vr_mn = zeros(size(x));
    vr_md = zeros(size(x));
end
