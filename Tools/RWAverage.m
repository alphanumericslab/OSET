function [mn, vr_mn, md, vr_md] = RWAverage(x)
%
% [mn, vr_mn, md, vr_md] = RWAverage(x)
% Robust weighted averaging
%
% inputs:
%   x: an (N x T) matrix containing N ensembles of a noisy event-related signal of length T
%
% output:
%   mn: the robust weighted average over the N rows of x
%   vr_mn: the variance of the average beat across the N rows of x
%   md: the robust weighted median over the N rows of x
%   vr_md: the variance of the median beat across the N rows of x
%
% Reference:
%   J.M. Leski. Robust weighted averaging of biomedical signals. IEEE Trans. Biomed. Eng., 49(8):796-804, 2002.
%
% Updates:
%   February 2019: return the average beat variance
%   November 2021: return the median beat and its variance
%
% Copyright Reza Sameni, Nov 2021
% The Open-Source Electrophysiological Toolbox
% (https://github.com/alphanumericslab/OSET)


if(size(x,1) > 1)
    % average beat (mean)
    mn0 = mean(x, 1);
    noise0 = x - mn0(ones(size(x, 1), 1), :);
    vr = var(noise0, [], 2);
    sm = sum(1./vr);
    weight = 1./(vr*sm);
    mn = weight' * x;
    noise = x - mn(ones(size(x, 1), 1), :);
    vr_mn = var(noise, [], 1);

    % average beat (median)
    md0 = median(x, 1);
    noise0 = x - md0(ones(size(x, 1), 1), :);
    vr = var(noise0, [], 2);
    sm = sum(1./vr);
    weight = 1./(vr*sm);
    md = weight' * x;
    noise = x - md(ones(size(x, 1), 1), :);
    vr_md = var(noise, [], 1);
else
    mn = x;
    md = x;
    vr_mn = zeros(size(x));
    vr_md = zeros(size(x));
end

