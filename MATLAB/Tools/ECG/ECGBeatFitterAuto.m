function [params model er indexes method] = ECGBeatFitterAuto(ECGmean,meanphase,approach,varargin),
%
% [params model er indexes method] = ECGBeatFitterAuto(ECGmean,meanphase,approach,stopth,divisions,maxit,energyth,wlen)
% Automatic fitting of sum of Gaussian kernels over a given ECG
%
% inputs:
% ECGmean: the input ECG beat
% meanphase: the phase signal of the ECG beat
% approach: 'fixed number' or 'min error'
% stopth: stopping threshold error percentage
% divisions: number of divisions for random point approach
% maxit: maximum number of iterations for randomized methods
% energyth: the energy threshold used for detecting peaks (default: 0.001)
% wlen: peak search window length (default: 0.01*SignalLength)
%
% output:
% params: the optimized parameters in the form of [alphai(1:N) bi(1:N) tetai(1:N)]
% model: the generated model (same length as ECGmean)
% er: the percentage of model error
% indexes: the locations of the Gaussian functions (before optimization)
% method: the approach used for finding the peaks
%
% Note:
%
% Open Source ECG Toolbox, version 2.1, September 2010
% Released under the GNU General Public License
% Copyright (C) 2010  Reza Sameni
% Shiraz University, Shiraz, Iran
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% beat length
N = length(ECGmean);

% idle segment at the beginning of the beat
L0 = ceil(N/10);

stopth = varargin{1};
num = varargin{2};
maxitr = varargin{3};
th = varargin{4};

if(nargin < 7 || isempty(varargin{4}))
    energyth = 0.001;
else
    energyth = varargin{4};
end

if(nargin < 8 || isempty(varargin{5}))
    wlen = ceil(N/100);
else
    wlen = varargin{5};
end

%//////////////////////////////////////////////////////////////////////////
% preprocessing
base = median(ECGmean(1:L0));
x = ECGmean - base;
y = BaseLine1(x,5,'md');%LPFilter(d,25/N);
% y = wden(x,'rigrsure','s','mln',6,'coif5');
y = LPFilter(y,45/N);

% % % figure;
% % % hold on;
% % % plot(x);
% % % plot(y,'r');
% % % grid

if(strcmp(approach,'min error'))

    %//////////////////////////////////////////////////////////////////////////
    % zero crossing detection
    zc = [0 y(1:end-1).*y(2:end)];
    I = find(zc<=0);
    E = sum(x.^2);
    II = [];
    m0 = find(I>L0,1);
    m = m0;
    for i = m0:length(I),
        ind = I(m):I(i);
        if(sum(x(ind).^2)>th*E)
            if(isempty(II))
                II = [I(m) I(i)];
            else
                II = [II I(i)];
            end
            m = m + 1;
        end
    end

    %//////////////////////////////////////////////////////////////////////////
    % local peak detection

    pk = zeros(1,N);
    for i = L0:N-L0,
        ind = max(i-wlen,1):min(i+wlen,N);
        if(y(i)==max(y(ind)) && y(i)>0 && std(y(ind))~=0) % local max with positive amplitudes
            pk(i) = 1;
        elseif(y(i)==min(y(ind)) && y(i)<0 && std(y(ind))~=0) % local min with negative amplitudes
            pk(i) = 1;
        end
    end
    J = find(pk);

    E = sum(y.^2);
    JJ = [];
    m0 = find(J>L0,1);
    m = m0;
    for i = m0:length(J),
        ind = J(m):J(i);
        if(sum(y(ind).^2)>th*E)
            if(isempty(JJ))
                JJ = [J(m) J(i)];
            else
                JJ = [JJ J(i)];
            end
            m = m + 1;
        end
    end

    %//////////////////////////////////////////////////////////////////////////
    % merge the results
    K = sort([I J]);
    KK = sort([II JJ]);

    % uniform points
    UU = L0:round(N/num):N-L0;

    %//////////////////////////////////////////////////////////////////////////
    % start optimizing
    done = logical(0);

    % peaks alone
    if(~done)
        indx = JJ;
        [params model er] = ECGOptimizeModel(y, indx, meanphase);
        indexes = indx;
        if(er<stopth)
            method = 0;
            done = 1;
        end
    end

    % peaks plus some random deviations
    if(~done)
        for i = 1:maxitr,
            indx = JJ + round((rand(size(JJ))-.5)*N/20);
            indx = sort(indx);
            indx = indx(indx>=1);
            indx = indx(indx<=N);
            [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
            if(er_<er)
                indexes = indx;
                params = params_;
                model = model_;
                er = er_;
                if(er<stopth)
                    method = 7;
                    done = 1;
                    break;
                end
            end
        end
    end

    % zero crossings alone
    if(~done)
        indx = II;
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 1;
                done = 1;
            end
        end
    end

    % zero-crossings and peaks together
    if(~done)
        indx = KK;
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 2;
                done = 1;
            end
        end
    end

    % uniform points
    if(~done)
        indx = UU;
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 3;
                done = 1;
            end
        end
    end

    % zero-crossings and their mid-points together
    if(~done)
        indx = sort([II round((II(1:end-1)+II(2:end))/2)]);
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 4;
                done = 1;
            end
        end
    end

    % peaks and their mid-points together
    if(~done)
        indx = sort([JJ round((JJ(1:end-1)+JJ(2:end))/2)]);
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 5;
                done = 1;
            end
        end
    end

    % zero-crossings, peaks, and their mid-points together
    if(~done)
        indx = sort([KK round((KK(1:end-1)+KK(2:end))/2)]);
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 6;
                done = 1;
            end
        end
    end

    % peaks and their mid-points plus some random deviations
    if(~done)
        for i = 1:maxitr,
            indx = [JJ round((JJ(1:end-1)+JJ(2:end))/2)];
            indx = indx + round((rand(size(indx))-.5)*N/10);
            indx = sort(indx);
            indx = indx(indx>=1);
            indx = indx(indx<=N);
            [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
            if(er_<er)
                indexes = indx;
                params = params_;
                model = model_;
                er = er_;
                if(er<stopth)
                    method = 8;
                    done = 1;
                    break;
                end
            end
        end
    end

    % zero-crossings and peaks plus some random deviations
    if(~done)
        for i = 1:maxitr,
            indx = KK + round((rand(size(KK))-.5)*N/10);
            indx = sort(indx);
            indx = indx(indx>=1);
            indx = indx(indx<=N);
            [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
            if(er_<er)
                indexes = indx;
                params = params_;
                model = model_;
                er = er_;
                if(er<stopth)
                    method = 9;
                    done = 1;
                    break;
                end
            end
        end
    end

    % raw zero-crossings alone
    if(~done)
        indx = I;
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 10;
                done = 1;
            end
        end
    end

    % raw peaks alone
    if(~done)
        indx = J;
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 11;
                done = 1;
            end
        end
    end

    % random points plus some random deviations
    if(~done)
        for i = 1:maxitr,
            indx = UU + round((rand(size(UU))-.5)*N/10);
            indx = sort(indx);
            indx = indx(indx>=1);
            indx = indx(indx<=N);
            [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
            if(er_<er)
                indexes = indx;
                params = params_;
                model = model_;
                er = er_;
                if(er<stopth)
                    method = 12;
                    done = 1;
                    break;
                end
            end
        end
    end

    if(~done)
        disp('Info: The specified requirements have not been met with fixed error model; continuing with random point approach.');
        approach = 'fixed number';
    end
elseif(strcmp(approach,'fixed number'))
    %//////////////////////////////////////////////////////////////////////////
    % uniform points
    UU = L0:round((N-2*L0)/num):N-L0;
    er = 100;
    for i = 1:maxitr,
        indx = UU + round((rand(size(UU))-.5)*N/10);
        indx = sort(indx);
        indx = indx(indx>=1);
        indx = indx(indx<=N);
        [params_ model_ er_] = ECGOptimizeModel(y, indx, meanphase);
        if(er_<er)
            indexes = indx;
            params = params_;
            model = model_;
            er = er_;
            if(er<stopth)
                method = 13;
                done = 1;
                break;
            end
        end
    end
end
if(~done)
    method = 14;
end
disp(['Final model fitting error percentage: ',num2str(er)]);
%//////////////////////////////////////////////////////////////////////////

% plot the results
n = 1:length(x);
figure;
plot(n,x);
hold on;
plot(n,y,'r');
plot(n,model,'g');
plot(n(indexes),x(indexes),'ro','linewidth',2);
grid

% % % grid
% % % disp(['Error = ',num2str(er),'%']);

% % % plot(d,'g');
% % % plot(dd,'m');
% % % plot(ddd,'k');
% % % plot(n(I),y(I),'ro','linewidth',3);
% % % plot(n(II),y(II),'go','linewidth',3);
% % % plot(n(J),y(J),'ro','linewidth',3);
% % % plot(n(JJ),y(JJ),'go','linewidth',3);
% % % plot(n(K),y(K),'ro','linewidth',1);
% % % plot(n(KK),y(KK),'go','linewidth',3);
% % % plot(n(UU),y(UU),'go','linewidth',3);
% % % grid
