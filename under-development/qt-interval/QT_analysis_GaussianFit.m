function [QT_ML, QT_BYS] = QT_analysis_GaussianFit(ecg, fs, varargin)

if nargin > 2
    soi1 = varargin{1};
else
    %% set the input values
    soi1.q = [-0.080; -0.020];
    soi1.t = [.1; .5];
end

if nargin > 3
    varNoise = varargin{2};
else
    varNoise = 0;
end

if nargin > 4
    plot_results = varargin{3};
else
    plot_results = false;
end

%% prepare data
L = size(ecg,1); % ecg length
tst = (0 : L-1)/fs; % time stamp
w1 = .75; % window length of the first median filter used in base line removing
w2 = .9; % window length of the second median filter used in base line removing
ecg = ecg - (BaseLine1(BaseLine1(ecg', round(w1*fs), 'md'), round(w2*fs), 'mn'))'; % baseline removal

for chnl = 1 : size(ecg, 2)

    rPeaks = find(PeakDetection20(ecg(:,chnl),70/60/fs,.6));


    [~, Q,~, ~, T] = ecgWavesSoI_RR(ecg(:,chnl), rPeaks, 1);
    soi2.q = Q(:,[1,3])/fs;
    soi2.t = T(:,[1,3])/fs;

    [~, Q,~, ~, T] = ecgWavesSoI_k(ecg(:,chnl), rPeaks, 1);
    soi3.q = Q(:,[1,3])/fs;
    soi3.t = T(:,[1,3])/fs;
    %% initial params
    % p0.q=[0; .02;  -.05];
    % p0.t=[0; .05; .25];
    p0=[];

    %% bounds
    lb=[]; ub=[];

    %% beta
    beta = 3;

    %% optimization setting
    options = struct('SpecifyObjectiveGradient',true);


    %% qt interval estimation
    qtInt1 = qtIntGausFit_v1(ecg(:,chnl), fs, rPeaks, 3, 3);
    qtInt2 = qtIntGausFit_v2(ecg(:,chnl), fs, rPeaks, 3, 3);
    qtInt3 = qtIntGausFit_v3(ecg(:,chnl), fs, rPeaks, 3, 3);

    %% ml framework
    soi = soi1;
    [mlGaussParams, rPeaks, soi,  waveParams, qtInt]=qtParamsGausFit_cl(ecg(:,chnl), fs, rPeaks, p0, soi, lb, ub, beta, options);

    if plot_results
        % polt the evaluated gaussians on the signal
        figure; p1=plot(tst,ecg(:,chnl)); % plot the channel
        hold on

        for j=1:length(rPeaks)
            tq = soi.q(j,1):1/fs:soi.q(j,2);
            tt = soi.t(j,1):1/fs:soi.t(j,2);
            p2 = plot(tq+rPeaks(j)/fs,GausVal(tq,mlGaussParams.q(:,j)),'r-');
            p3 = plot(tt+rPeaks(j)/fs,GausVal(tt,mlGaussParams.t(:,j)),'r-');
            p4 = plot([waveParams.q(1,j) waveParams.t(2,j)]+rPeaks(j)/fs, ...
                ecg(floor(fs*[waveParams.q(1,j) waveParams.t(2,j)])+rPeaks(j),chnl),'c*');
        end
        legend([p1 p2 p4],'ecg', 'Gaussians', 'q/t onset/offset')
        title 'ML framework'
    end

    %% Bys framework
    PrMu.q=mean(mlGaussParams.q(:,:),2);
    PrCov.q=cov(mlGaussParams.q(:,:)');

    PrMu.t=mean(mlGaussParams.t(:,:),2);
    PrCov.t=cov(mlGaussParams.t(:,:)');

    [bysGaussParams, rPeaks, soi, waveParams, qtInt]=qtParamsGausFit_cl(ecg(:,chnl), fs, rPeaks, p0, soi, lb, ub, beta, PrMu, PrCov, varNoise);

    if plot_results
        % polt the evaluated gaussians on the signal
        figure; p1=plot(tst,ecg(:,chnl)); % plot the channel
        hold on

        for j=1:length(rPeaks)
            tq=soi.q(j,1):1/fs:soi.q(j,2);
            tt=soi.t(j,1):1/fs:soi.t(j,2);
            p2=plot(tq+rPeaks(j)/fs,GausVal(tq,bysGaussParams.q(:,j)),'r-');
            p3=plot(tt+rPeaks(j)/fs,GausVal(tt,bysGaussParams.t(:,j)),'r-');
            p4=plot([waveParams.q(1,j) waveParams.t(2,j)]+rPeaks(j)/fs, ...
                ecg(floor(fs*[waveParams.q(1,j) waveParams.t(2,j)])+rPeaks(j),chnl),'c*');
        end
        legend([p1 p2 p4],'ecg', 'Gaussians', 'q/t onset/offset')
        title 'BYS framework'
    end
end

    %{
function [QT_ML, QT_BYS] = QT_analysis_GaussianFit(ecg, fs, varargin)

if nargin > 2
    soi = varargin{1};
else
    %% set the input values
    soi.q=[-0.045; -0.015];
    soi.t=[.1; .5];
end

if nargin > 3
    varNoise = varargin{2};
else
    varNoise = 0;
end

if nargin > 4
    plot_results = varargin{3};
else
    plot_results = false;
end

L = size(ecg, 1); % ecg length
tst = (0 : L-1) / fs; % time stamp

%% ml framework
[GaussParams_ML, rPeaks_ML, soi_ML,  waveParams_ML, QT_ML] = qtParamsGausFit(ecg, fs, [],[], soi);

for chnl = 1 : size(ecg, 2)
    PrMu.q(:, chnl) = mean(GaussParams_ML.q(:,:,chnl),2);
    PrCov.q(:, :, chnl) = cov(GaussParams_ML.q(:,:,chnl)');

    PrMu.t(:, chnl) = mean(GaussParams_ML.t(:,:,chnl),2);
    PrCov.t(:, :, chnl) = cov(GaussParams_ML.t(:,:,chnl)');
end

%% Bayesian framework
[GaussParams_BYS, rPeaks_BYS, soi_BYS, waveParams_BYS, QT_BYS] = qtParamsGausFit(ecg, fs, [], [], soi_ML, [], [], [], PrMu, PrCov, varNoise);


if plot_results
    for chnl = 1 : size(ecg, 2)

        % plot the evaluated gaussians on the signal
        figure;
        p1 = plot(tst,ecg(:,chnl), 'b'); % plot the channel
        hold on
        tq = soi_ML.q(1):1/fs:soi_ML.q(2);
        tt = soi_ML.t(1):1/fs:soi_ML.t(2);
        for j = 1:length(rPeaks_ML)
            p2 = plot(tq + rPeaks_ML(j)/fs, GausVal(tq, GaussParams_ML.q(:,j,chnl)),'k-');
            p3 = plot(tt + rPeaks_ML(j)/fs, GausVal(tt, GaussParams_ML.t(:,j,chnl)),'r-');
            QT_indexes_ML = floor(fs*[waveParams_ML.q(1, j, chnl) waveParams_ML.t(2, j, chnl)] + rPeaks_ML(j));
            QT_indexes_ML(QT_indexes_ML < 1) = nan;
            QT_indexes_ML(QT_indexes_ML > length(ecg(:,chnl))) = nan;
            p4 = plot(QT_indexes_ML/fs, ecg(QT_indexes_ML, chnl),'c*');
        end
        legend([p1, p2 p3, p4],'ecg', 'Q-wave fit', 'T-wave fit', 'q/t onset/offset')
        title 'ML framework'

        figure;
        p1 = plot(tst,ecg(:,chnl), 'b'); % plot the channel
        hold on
        tq = soi_BYS.q(1):1/fs:soi_BYS.q(2);
        tt = soi_BYS.t(1):1/fs:soi_BYS.t(2);
        for j = 1:length(rPeaks_BYS)
            p2 = plot(tq + rPeaks_BYS(j)/fs, GausVal(tq,GaussParams_BYS.q(:,j)), 'k-');
            p3 = plot(tt + rPeaks_BYS(j)/fs, GausVal(tt,GaussParams_BYS.t(:,j)), 'r-');
            QT_indexes_BYS = floor(fs*[waveParams_BYS.q(1, j, chnl) waveParams_BYS.t(2, j, chnl)] + rPeaks_BYS(j));
            QT_indexes_BYS(QT_indexes_BYS < 1) = nan;
            QT_indexes_BYS(QT_indexes_BYS > length(ecg(:,chnl))) = nan;
            p4 = plot(QT_indexes_BYS/fs, ecg(QT_indexes_BYS, chnl),'c*');
        end
        legend([p1, p2 p3, p4],'ecg', 'Q-wave fit', 'T-wave fit', 'q/t onset/offset')
        title 'BYS framework'
    end
end
    %}
