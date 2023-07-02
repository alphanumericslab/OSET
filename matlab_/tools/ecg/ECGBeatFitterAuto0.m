function [prms,mdl,error,approach] = ECGBeatFitterAuto(ECGmean,meanphase),
% clc
% clear all
% close all;
%
% load('SampleECG2.mat'); data = data(1:15000,8)';
%
% fs = 1000;
% t = [0:length(data)-1]/fs;
%
% f = 1;                                          % approximate R-peak frequency
%
% bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
% x = data-bsline;
% peaks = PeakDetection(x,f/fs);                  % peak detection
%
% [phase phasepos] = PhaseCalculation(peaks);     % phase calculation
%
% teta = 0;                                       % desired phase shift
% pphase = PhaseShifting(phase,teta);             % phase shifting
%
% bins = 500;                                     % number of phase bins
% [ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1);

bins = length(ECGmean);
sm = BaseLine1(ECGmean,5,'mn');

localpeaks = zeros(6,bins);
% localpeaks = zeros(4,bins);
prd = [1 sm(1:end-1).*sm(2:end)];
prd(1) = -1;
prd(end) = -1;
zerocrs = logical(prd<=0);
II = find(zerocrs);
S = sum(sm.^2);
for j = 1:length(II)-1,
    indj = II(j):II(j+1);
    Sj = sum(sm(indj).^2);
    if(Sj>=(.0002*S))
        winjlen = length(indj);
        [mmx maxj] = max(abs(sm(indj)));
        mmd = median(sm(indj));
        mmdj = find(sm(indj)==mmd);

        % select the average of the two ending points of the segment
        localpeaks(1,round(mean([indj(1),indj(end)]))) = 1;

        % select the peak point of the segment
        localpeaks(2,maxj+indj(1)-1) = 1;

        if(abs(maxj-winjlen/2)>floor(winjlen/10)),
            localpeaks(3,indj(1)-1 + ceil(winjlen/3)) = 1;
            localpeaks(3,indj(end) - ceil(winjlen/3)) = 1;
        else
            localpeaks(3,indj(1)-1 + maxj) = 1;
        end

        if(abs(maxj-winjlen/2)>floor(winjlen/10)),
            localpeaks(4,indj(1)-1 + maxj) = 1;
            localpeaks(4,round(mean([indj(1),indj(end)]))) = 1;
        else
            localpeaks(4,indj(1)-1 + maxj) = 1;
        end

        if(abs(maxj-winjlen/2)>floor(winjlen/10)),
            localpeaks(5,indj(1)-1 + ceil(winjlen/4)) = 1;
            localpeaks(5,indj(1)-1 + ceil(winjlen/2)) = 1;
            localpeaks(5,indj(1)-1 + ceil(3*winjlen/4)) = 1;
        else
            localpeaks(5,indj(1)-1 + maxj) = 1;
        end

        if(abs(maxj-winjlen/2)>floor(winjlen/10)),
            localpeaks(6,indj(1)-1 + ceil(winjlen/4)) = 1;
            localpeaks(6,indj(1)-1 + ceil(winjlen/2)) = 1;
            localpeaks(5,indj(1)-1 + ceil(3*winjlen/4)) = 1;
        else
            localpeaks(6,max(indj(1)-1 + maxj-ceil(winjlen/4),indj(1))) = 1;
            localpeaks(6,min(indj(1)-1 + maxj+ceil(winjlen/4),indj(end))) = 1;
        end
    end
end

Model = zeros(size(localpeaks,1),bins);
er = zeros(size(localpeaks,1),1);
for i = 1:size(localpeaks,1),
    I = find(localpeaks(i,:));
    P = length(I);

    tetai = meanphase(I(1:P));
    alphai = 1*sm(I(1:P));
    bi = .01*ones(size(alphai));

    options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100);
    InitParams = [alphai bi tetai];

    OptParams{i} = nlinfit(meanphase,sm,@ECGModel,InitParams,options);
    % OptParams = lsqnonlin(@(InitParams) ECGModelError(InitParams,ECGmn,Phasemn,0),InitParams,InitParams-2,InitParams+2,options);
    % Model0 = ECGModelError(InitParams,sm,meanphase,1);
    Model(i,:) = ECGModelError(OptParams{i},sm,meanphase,1);
    er(i) = 100*mean((sm-Model(i,:)).^2)/mean(sm.^2);
    %plot(meanphase,Model(i,:),sty(i));
end

[Y,II] = min(er);
prms = OptParams{II};
mdl = Model(II,:);
error = er(II);
approach = II;
