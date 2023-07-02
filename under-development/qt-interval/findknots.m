function k=findknots(x,nknots, varargin)
% 
% Syntax:
% k=ecgknots(ecg, nknots)
% k=ecgknots(ecg, nknots, dist)
% 
% Description:
% Finding the suitable knots along the ecg signal for spline or polynomial
% fitting. At the first level (depth=0), the first and last samples are
% set as initial knots, and for each next levels, the point with maximum 
% distance from the reference line (the line passing through the knots) are 
% considered as a new knot.
% 
% Inputs:
%   ecg: A beat of ecg signal.
%   nknots: number of desired knots
%   dist: distance measure; accepts 'Euclidean' or 'Value'. 
%   
% 
% Output:
%   k: Vector of the found knots
% 
% 
% Davood Fattahi, 10/02/2020;
% fattahi.d@gmail.com
% 
% update 24/02/2022: 
%   1) optional input of 'nknots' is added. 
%   3) unnecessary input parameters (depth and threshold) are removed.
%   2) bug of adjacent knots (one sample distanced) is resolved.
%   3) line slope formula is modified.
%   4) optional inputs of value-distance and Euclidean distance are added. -Davood Fattahi
% 

if isempty(varargin)
    Ecd = true; Vald = false;
elseif strcmp(varargin{1},'Euclidean')
    Ecd = true; Vald = false;
elseif strcmp(varargin{1},'Value')
    Ecd = false; Vald = true;
end
if nknots<3
    error('number of knots should be greather than 2')
end

x=x(:)';
k=nan(nknots-1,2);
k(1,1)=1;
k(1,2)=size(x(:),1);
for j=2:nknots-1
    M = nan(j-1,1); I = nan(j-1,1);
    for i=1:j-1
        tt{i}=k(i,1):k(i,2);
        s=x(tt{i});
        if Ecd
            m = (s(end)-s(1)) / (tt{i}(end)-tt{i}(1));
            [M(i),I(i)] = max(abs((-m*tt{i}+s)+m*tt{i}(1)-s(1))./sqrt(m^2+1));
        elseif Vald
            line = s(1) + (0:length(s)-1).*((s(end)-s(1))/(length(s)-1));
            [M(i),I(i)] = max(abs(s-line)); %%% value distance
        end
        
   
    end
    [~,II] = max(M);
    k(j,1) = tt{II}(I(II));
    k(j,2) = k(II,2);
    k(II,2) = tt{II}(I(II));
end

k=unique(k(:));

return

