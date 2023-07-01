%
% Test program for mean ECG extraction.
%
% By:
% Reza Sameni, September 2006
% LIS-INPG, Grenoble, France - Sharif University of Technology, Tehran, Iran
% reza.sameni@gmail.com

clc
clear
close all;

load('SampleECG2.mat'); data = data(1:15000,8)';

fs = 1000;
t = [0:length(data)-1]/fs;

f = 1;                                          % approximate R-peak frequency

bsline = LPFilter(data,.7/fs);                  % baseline wander removal (may be replaced by other approaches)
x = data-bsline;
peaks = PeakDetection(x,f/fs);                  % peak detection

[phase phasepos] = PhaseCalculation(peaks);     % phase calculation

teta = 0;                                       % desired phase shift
pphase = PhaseShifting(phase,teta);             % phase shifting

bins = 500;                                     % number of phase bins
[ECGmean,ECGsd,meanphase] = MeanECGExtraction(x,pphase,bins,1);

sm = BaseLine1(ECGmean,5,'mn');
% df1 = [sm(1) diff(sm)];
% df1 = BaseLine1(df1,10,'mn');
% df2 = [df1(1) diff(df1)];
% sm = cumsum(ECGmean);

localpeaks = zeros(6,bins);
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
            localpeaks(3,indj(1)-1 + floor(winjlen/3)) = 1;
            localpeaks(3,indj(end) - floor(winjlen/3)) = 1;
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
            localpeaks(5,indj(1)-1 + floor(winjlen/4)) = 1;
            localpeaks(5,indj(1)-1 + floor(winjlen/2)) = 1;
            localpeaks(5,indj(1)-1 + floor(3*winjlen/4)) = 1;
        else
            localpeaks(5,indj(1)-1 + maxj) = 1;
        end

        if(abs(maxj-winjlen/2)>floor(winjlen/10)),
            localpeaks(6,indj(1)-1 + floor(winjlen/4)) = 1;
            localpeaks(6,indj(1)-1 + floor(winjlen/2)) = 1;
            localpeaks(5,indj(1)-1 + floor(3*winjlen/4)) = 1;
        else
            localpeaks(6,max(indj(1)-1 + maxj-floor(winjlen/4),indj(1))) = 1;
            localpeaks(6,min(indj(1)-1 + maxj+floor(winjlen/4),indj(end))) = 1;
        end
        % % %         %if(Sj>=(.1*S) || abs(maxj-winjlen/2)>floor(winjlen/10)),
        % % %         if(abs(maxj-winjlen/2)>floor(winjlen/10)),
        % % %             %         if(abs(mmdj-winjlen/2)>floor(winjlen/5000)),
        % % %             %if(abs(mmdj+indj(1)-1 - round(mean([indj(1),indj(end)])))>=.05*winjlen)
        % % %             %         if(abs(mmdj - maxj)>=.00001*winjlen)
        % % %
        % % %             %             peaks4(indj(1)-1 + floor(winjlen/4)) = 1;
        % % %             %             peaks4(indj(1)-1 + floor(winjlen/2)) = 1;
        % % %             %             peaks4(indj(1)-1 + floor(3*winjlen/4)) = 1;
        % % %
        % % %             peaks4(indj(1)-1 + floor(winjlen/3)) = 1;
        % % %             peaks4(indj(end) - floor(winjlen/3)) = 1;
        % % %             %                         peaks4(indj(1)-1 + maxj) = 1;
        % % %             %                         peaks4(round(mean([indj(1),indj(end)]))) = 1;
        % % %         else
        % % %             peaks4(indj(1)-1 + maxj) = 1;
        % % %         end

    end
end

sty = ['b','r','g','m','c','y'];
figure;
plot(meanphase,sm,'k','linewidth',2);
hold on;
xlabel('phase (rads.)');
grid;

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
    plot(meanphase,Model(i,:),sty(i));
end

% % % % % figure;
% % % % % plot(t,data*6/max(data),'b');
% % % % % hold on
% % % % % plot(t,peaks*2,'ro');
% % % % % plot(t,phase,'c--','linewidth',1);
% % % % % plot(t,pphase,'g','linewidth',1);
% % % % % grid;
% % % % % xlabel('time (sec.)');
% % % % % legend('Scaled ECG','ECG Peaks','phase','shifted phase');
% % % % % 
% % % % % figure;
% % % % % %errorbar(meanphase,ECGmean,ECGsd/2);
% % % % % hold on;
% % % % % plot(meanphase,ECGmean,'r');
% % % % % % plot(meanphase,peaks2/2,'y','linewidth',2);
% % % % % % plot(meanphase,peaks3/2,'c','linewidth',2);
% % % % % 
% % % % % plot(meanphase,peaks4/2,'k','linewidth',1);
% % % % % 
% % % % % plot(meanphase,sm,'b');
% % % % % % plot(meanphase,Model0,'m');
% % % % % plot(meanphase,Model,'g');
% % % % % 
% % % % % % plot(meanphase,df1*10,'c','linewidth',2);
% % % % % % plot(meanphase,df2*30,'y','linewidth',2);
% % % % % % plot(meanphase,zerocrs,'bo');
% % % % % 
% % % % % % legend('SD Bar','Mean ECG');
