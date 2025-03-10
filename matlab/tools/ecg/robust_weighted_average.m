function [mn, vr_mn, varargout] = robust_weighted_average(x)
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
%   md: The robust weighted median over the N rows of x (optional)
%   vr_md: The variance of the median beat across the N rows of x (optional)
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
%   2023: Added varargout to speed up when median estimates are not required
%   2025: Handling NAN values in data
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

num_beats = size(x, 1);
if num_beats > 1
    % Average beat (mean)
    mn0 = mean(x, 1,'omitnan');
    noise0 = x - mn0;
    vr = var(noise0, [], 2,'omitnan');
    sm = sum(1 ./ vr);
    weight = 1 ./ (vr * sm);
    if any(isnan(x(:)))
        weight_mat = repmat(weight,1,size(x, 2));
        weight_mat(isnan(x))=0;
        ind_nan = sum(weight_mat,1)==0;
        weight_mat = weight_mat./sum(weight_mat,1);
        mn = sum(weight_mat.* x,1,'omitnan');
        mn(ind_nan) = nan;
        noise = x - mn;
        vr_mn = var(noise, [], 1,'omitnan');
    else
        mn = weight'* x;
        noise = x - mn;
        vr_mn = var(noise, [], 1,'omitnan');
    end

    if nargout > 2
        % Average beat (median)
        md0 = median(x, 1,'omitnan');
        noise0 = x - md0;
        vr = var(noise0, [], 2,'omitnan');
        sm = sum(1 ./ vr);
        weight = 1 ./ (vr * sm);
        if any(isnan(x(:)))
            weight_mat = repmat(weight,1,size(x, 2));
            weight_mat(isnan(x))=0;
            ind_nan = sum(weight_mat,1)==0;
            weight_mat = weight_mat./sum(weight_mat,1);
            md = sum(weight_mat.* x,1,'omitnan');
            md(ind_nan) = nan;

        else
            md = weight'* x;
        end

        varargout{1} = md;
        if nargout > 3
            noise = x - md;
            vr_md = var(noise, [], 1,'omitnan');
            varargout{2} = vr_md;
        end
    end

else
    mn = x;
    vr_mn = zeros(size(x));
    if nargout > 2
        md = x;
        varargout{1} = md;
    end

    if nargout > 3
        vr_md = zeros(size(x));
        varargout{2} = vr_md;
    end
end

end
