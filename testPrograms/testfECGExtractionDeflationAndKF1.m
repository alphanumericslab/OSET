% double KF for fetal denoising

clear;
close all;

load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; time = FOETAL_ECG(:,1)'; clear FOETAL_ECG; fs = 250;

L1 = size(data,1);
L2 = size(data,2);
bins = fs;

% R-peak detection
ref = data(1,:);
b = LPFilter(ref,1.5/fs);
ref = ref - b;
f = 1.1;
peaks = PeakDetection(ref,f/fs);

% ECG phase calculation
[phase, phase2] = PhaseCalculation(peaks);

% PM time calculation
J = find(peaks);
n1 = diff(J);
prd = round(mean(n1));
wlen = max(n1)-min(n1);

T1 = zeros(L2-prd-wlen,1);
NN = length(T1);
for t = 1:NN
    df = abs(phase(t) - phase(t+prd-wlen:t+prd+wlen));
    [Y,I] = min(df);
    T1(t) = t + prd + I - wlen -1;
end
T0 = 1:NN;

% Main iteration
Itr = 5;
wholedata = zeros(L1,L2,Itr+1);
wholedata(:,:,1) = data;
dat = data;

% BW removal
bl = LPFilter(dat,.7/fs);
dat = dat - bl;

W0 = jadeR(dat);
s0 = W0*dat;
s0(1:6,:) = 0;
dat0 = pinv(W0)*s0;
% % % PlotECG(dat0,8,'b');

%//////////////////////////////////////////////////////////////////////////
for i = 1:Itr
    % periodic component analysis stage
    A = dat(:,T0)*dat(:,T1)';
    B = dat(:,T0)*dat(:,T0)';

    A = (A+A')/2;
    B = (B+B')/2;

    [V,D] = eig(A,B,'chol');

    d = diag(D);
    [Y,I] = sort(d);
    I = I(end:-1:1);

    W = V(:,I)';
    dat = W*dat;

    % remove the ECG from the most peridic component
    xx = dat(1,:);
    %//////////////////////////////////////////////////////////////////////////
    bsline = LPFilter(xx,.7/fs);
    xxx = xx - bsline;
    [ECGmean,ECGsd,meanphase] = MeanECGExtraction(xxx,phase,bins,1); % mean ECG extraction

    ECGBeatFitter(ECGmean,ECGsd,meanphase);           % ECG beat fitter GUI

    %//////////////////////////////////////////////////////////////////////////
    N = length(OptimumParams)/3;% number of Gaussian kernels
    JJ = find(peaks);
    fm = fs./diff(JJ);          % heart-rate
    w = mean(2*pi*fm);          % average heart-rate in rads.
    wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.

    y = [phase ; xxx];

    X0 = [-pi 0]';
    P0 = [(2*pi)^2 0 ;0 (10*max(abs(xxx))).^2];
    Q = diag( [ (.1*OptimumParams(1:N)).^2 (.05*ones(1,N)).^2 (.05*ones(1,N)).^2 (wsd)^2 , (.03*mean(ECGsd(1:round(length(ECGsd)/10))))^2] );
    R = [(w/fs).^2/12 0 ;0 (mean(ECGsd(1:round(length(ECGsd)/10)))).^2];
    Wmean = [OptimumParams w 0]';
    Vmean = [0 0]';
    Inits = [OptimumParams w fs];

    InovWlen = ceil(.5*fs);     % innovations monitoring window length
    tau = [];                   % Kalman filter forgetting time. tau=[] for no forgetting factor
    gamma = 1;                  % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
    RadaptWlen = ceil(fs/2);    % window length for observation covariance adaptation

    %//////////////////////////////////////////////////////////////////////////
    [Xekf,Phat,Xeks,PSmoothed,ak] = EKSmoother(y,X0,P0,Q,R,Wmean,Vmean,Inits,InovWlen,tau,gamma,RadaptWlen,0);

    % returning back the residue
    dat(1,:) = xx - Xeks(2,:);

    dat = pinv(W)*dat;

    wholedata(:,:,i+1) = dat;
end
%//////////////////////////////////////////////////////////////////////////
% % % for i = 1:Itr+1,
% % %     %     if i == 1,
% % %     %         PlotECG(squeeze(wholedata(I,:,i)),8,'b');
% % %     %     else
% % %     PlotECG(squeeze(wholedata(:,:,i)),8,'b');
% % %     %     end
% % % end
% % %
% % % WW = jader(squeeze(wholedata(:,:,end)));
% % % s = WW*squeeze(wholedata(:,:,end));

% % % PM = zeros(Itr+1,L1);
% % % for itr = 1:Itr+1,
% % %     figure;
% % %     for i = 1:L1,
% % %         d = squeeze(wholedata(i,:,itr));
% % %         tmp = corrcoef(d(T0),d(T1));
% % %         PM(itr,i) = tmp(1,2);
% % %         subplot(L1,1,i);
% % %         plot(time,d,'k','linewidth',1);
% % %         grid;
% % %         set(gca,'Box','On','FontSize',16);
% % %         if (i<L1)
% % %             set(gca,'XTickLabel',[]);
% % %         end
% % %         %     ylabel(['IC_',num2str(i)],'FontSize',16);
% % %     end
% % %     xlabel('Time (s)','FontSize',16);
% % % end
% % %
% % % PM = round(PM*1000)/1000;
% % %
% % % PM

%//////////////////////////////////////////////////////////////////////////
maternalremoved = dat;

% R-peak detection
ref = dat(1,:);
b = LPFilter(ref,1.5/fs);
ref = ref - b;
f = 2;
fpeaks = PeakDetection(ref,f/fs);
% % % % fetal ECG phase calculation
% % % [fphase fphase2] = PhaseCalculation(fpeaks);
% % %
% % % [ECGmean,ECGsd,meanphase] = MeanECGExtraction(xxx,phase,bins,1); % mean ECG extraction

[s,W,A] = PiCA(maternalremoved,fpeaks);

PlotECG(squeeze(wholedata(:,:,end)),8,'b');

% PlotECG(s,8,'b');
s1 = s;
s1(3:end,:) = 0;
dat2 = A*s1;
% PlotECG(dat2,8,'b');

L1 = 8;
figure;
for i = 1:L1
    subplot(L1,1,i);
    hold on;
    plot(time,dat0(i,:),'b','linewidth',1);
    plot(time,dat2(i,:),'r','linewidth',1);
    grid;
    set(gca,'Box','On','FontSize',16);
    if (i<L1)
        set(gca,'XTickLabel',[]);
    end
    %     ylabel(['IC_',num2str(i)],'FontSize',16);
end
xlabel('Time (s)','FontSize',16);

% % % % %
% % % % % % ECG phase calculation
% % % % % [phase phase2] = PhaseCalculation(peaks);
% % % % %
% % % % % % PM time calculation
% % % % % J = find(peaks);
% % % % % n1 = diff(J);
% % % % % prd = round(mean(n1));
% % % % % wlen = max(n1)-min(n1);
% % % % %
% % % % % T1 = zeros(L2-prd-wlen,1);
% % % % % NN = length(T1);
% % % % % for t = 1:NN,
% % % % %     df = abs(phase(t) - phase(t+prd-wlen:t+prd+wlen));
% % % % %     [Y,I] = min(df);
% % % % %     T1(t) = t + prd + I - wlen -1;
% % % % % end
% % % % % T0 = 1:NN;
% % % % %
% % % % % % Main iteration
% % % % % Itr = 3;
% % % % % wholedata = zeros(L1,L2,Itr+1);
% % % % % wholedata(:,:,1) = dat;
% % % % %
% % % % % % BW removal
% % % % % bl = LPFilter(dat,.7/fs);
% % % % % dat = dat - bl;
% % % % %
% % % % % for i = 1:Itr,
% % % % %     % periodic component analysis stage
% % % % %     A = dat(:,T0)*dat(:,T1)';
% % % % %     B = dat(:,T0)*dat(:,T0)';
% % % % %
% % % % %     A = (A+A')/2;
% % % % %     B = (B+B')/2;
% % % % %
% % % % %     [V,D] = eig(A,B,'chol');
% % % % %
% % % % %     d = diag(D);
% % % % %     [Y,I] = sort(d);
% % % % %     I = I(end:-1:1);
% % % % %
% % % % %     W = V(:,I)';
% % % % %     dat = W*dat;
% % % % %
% % % % %     % remove the ECG from the most peridic component
% % % % %     xx = dat(1,:);
% % % % %     %//////////////////////////////////////////////////////////////////////////
% % % % %     bsline = LPFilter(xx,.7/fs);
% % % % %     xxx = xx - bsline;
% % % % %     [ECGmean,ECGsd,meanphase] = MeanECGExtraction(xxx,phase,bins,1); % mean ECG extraction
% % % % %
% % % % %     ECGBeatFitter(ECGmean,ECGsd,meanphase);           % ECG beat fitter GUI
% % % % %
% % % % %     %//////////////////////////////////////////////////////////////////////////
% % % % %     N = length(OptimumParams)/3;% number of Gaussian kernels
% % % % %     JJ = find(peaks);
% % % % %     fm = fs./diff(JJ);          % heart-rate
% % % % %     w = mean(2*pi*fm);          % average heart-rate in rads.
% % % % %     wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.
% % % % %
% % % % %     y = [phase ; xxx];
% % % % %
% % % % %     X0 = [-pi 0]';
% % % % %     P0 = [(2*pi)^2 0 ;0 (10*max(abs(xxx))).^2];
% % % % %     Q = diag( [ (.1*OptimumParams(1:N)).^2 (.05*ones(1,N)).^2 (.05*ones(1,N)).^2 (wsd)^2 , (.03*mean(ECGsd(1:round(length(ECGsd)/10))))^2] );
% % % % %     R = [(w/fs).^2/12 0 ;0 (mean(ECGsd(1:round(length(ECGsd)/10)))).^2];
% % % % %     Wmean = [OptimumParams w 0]';
% % % % %     Vmean = [0 0]';
% % % % %     Inits = [OptimumParams w fs];
% % % % %
% % % % %     InovWlen = ceil(.5*fs);     % innovations monitoring window length
% % % % %     tau = [];                   % Kalman filter forgetting time. tau=[] for no forgetting factor
% % % % %     gamma = 1;                  % observation covariance adaptation-rate. 0<gamma<1 and gamma=1 for no adaptation
% % % % %     RadaptWlen = ceil(fs/2);    % window length for observation covariance adaptation
% % % % %
% % % % %     %//////////////////////////////////////////////////////////////////////////
% % % % %     [Xekf,Phat,Xeks,PSmoothed,ak] = EKSmoother(y,X0,P0,Q,R,Wmean,Vmean,Inits,InovWlen,tau,gamma,RadaptWlen,0);
% % % % %
% % % % %     % returning back the residue
% % % % %     dat(1,:) = xx - Xeks(2,:);
% % % % %
% % % % %     dat = pinv(W)*dat;
% % % % %
% % % % %     wholedata(:,:,i+1) = maternalremoved - dat;
% % % % % end
% % % % %
% % % % % PM = zeros(Itr+1,L1);
% % % % % for itr = 1:Itr+1,
% % % % %     figure;
% % % % %     for i = 1:L1,
% % % % %         d = squeeze(wholedata(i,:,itr));
% % % % %         tmp = corrcoef(d(T0),d(T1));
% % % % %         PM(itr,i) = tmp(1,2);
% % % % %         subplot(L1,1,i);
% % % % %         plot(time,d,'k','linewidth',1);
% % % % %         grid;
% % % % %         set(gca,'Box','On','FontSize',16);
% % % % %         if (i<L1)
% % % % %             set(gca,'XTickLabel',[]);
% % % % %         end
% % % % %         %     ylabel(['IC_',num2str(i)],'FontSize',16);
% % % % %     end
% % % % %     xlabel('Time (s)','FontSize',16);
% % % % % end
% % % % % PM = round(PM*1000)/1000;
% % % % %
% % % % % PM
% % % % %
