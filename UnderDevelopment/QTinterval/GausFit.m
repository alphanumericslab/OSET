function [p]=GausFit(t, x, p0, lb, ub, varargin)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit sum of Gaussian functions on the data.
% 
% Syntax:
% [p, varargout]=GausFit(t, x, p0, lb, ub)
% [p, varargout]=GausFit(t, x, p0, lb, ub, options)
% [p, varargout]=GausFit(t, x, p0, lb, ub, PrMu, PrCov, PstVar)
% [p, varargout]=GausFit(t, x, p0, lb, ub, PrMu, PrCov, PstVar, options)
% 
% Inputs:
% t: time samples
% x: input data
% p0: parameters' initial values. The first, second, and third rows contain
%   respectively Gaussians' amplitude (a), Gaussians' width (b), and
%   Gaussians' center (c). So the i-th column contains the i-th Gaussian
%   parameters, by the order of [ai ; bi ; ci]. It also can be in the
%   vector form, like [a1 b1 c1 a2 b2 c2 ...]'.
% lb: lower bound for parameters.
% ub: upper bound for parameters.
% options: options structure for optimization problem. (Note: It should be in a
%   simple strutre format, NOT the 'optimoptions' format commonly used in matlab
%   optimization toolbox. However, the fields are same as optimoptions.)
% PrMu: mean of prior distribution.
% PrCov: covarinace matrix of prior information.
% PstVar: variance of the additive noise (scalar).
% 
% Output:
% p: the estimated parameters of Gaussians, with the same format as p0;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021  Davood Fattahi
% fattahi.d@gmail.com
% Developed as a part of implementation in https://doi.org/10.36227/techrxiv.15176139.v2
% 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

p0=p0(:); t=t(:); x=x(:); lb=lb(:); ub=ub(:);
NumGaus=length(p0)/3;
if round(NumGaus)~=NumGaus
    error('Wrong number of parameters in p0!')
end

if length(varargin)<=1 % ML
    ML=true; Bys=false;
elseif length(varargin)>=3 % Bys
    ML=false; Bys=true;
    PrMu=varargin{1}; PrMu=PrMu(:);
    PrCov=varargin{2};
    PstVar=varargin{3};
end
if nargin==5 || nargin==8
    options0=struct;
elseif nargin==6 || nargin==9
    options0=varargin{end};
end


%%% Gaussfit

if ML % ML
    fun=@(p)(MlFun(t,p,x));
elseif Bys % Bys
    fun=@(p)(BysFun(t,p,x,PstVar, PrCov, PrMu));
end
options = optimoptions(@lsqnonlin,'Display','off');
options=optionsmerge(options, options0);
p = lsqnonlin(fun,p0,lb,ub,options);

end

%% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%%% functions

function [F, J]=MlFun(t,p,x)
    F=GausVal(t,p)-x;
    if nargout>1
        J=GausGrad(t,p)';
    end
end

function [F, J]=BysFun(t,p,x,PstVar, PrCov, PrMu)
    F=([(GausVal(t,p)-x)./sqrt(PstVar) ; sqrtm(PrCov)\(p-PrMu)]);
    if nargout>1
        J=[GausGrad(t,p)'./sqrt(PstVar) ; inv(sqrtm(PrCov))];
    end
end

function options=optionsmerge(options1, options2)
%%% \\\\\\\\\\\\\\\\\\\\\\
%%% merging options2 into options1.
FN=fieldnames(options2);
for i=1:length(FN)
    eval(['options1.' FN{i} '=options2.' FN{i} ';'])
end
options=options1;
end
