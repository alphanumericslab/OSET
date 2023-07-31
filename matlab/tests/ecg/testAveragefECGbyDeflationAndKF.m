% test fetal ECG morphology extraction by the deflation algorithm and double KF for fetal denoising
% Copyright Reza Sameni, 2018
% Created 2008
% Modified June 2018
%
% Modified Oct 2020:
%   ECGBeatFitter replaced with new function call format

clear;
close all;

load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; time = FOETAL_ECG(:,1)'; clear FOETAL_ECG; fs = 250;

L1 = size(data,1);
L2 = size(data,2);
bins = fs/2;

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


% BW removal
data = data - LPFilter(data,.7/fs);
% for i = 1:size(data,1),
%     data(i,:) = data(i,:) - LPFilter(Median(data(i,:),size(data,2),round(.15*fs),round(.3*fs))',10/fs);
% end

% Main iteration
Itr = 4;
dat = data;
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
    
    OptimumParams = ECGBeatFitter(ECGmean,ECGsd,meanphase);           % ECG beat fitter GUI
    
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
    
    data2 = dat;
end

fetref = data2(1,:);
fpeaks = PeakDetection(fetref,2.2/fs);
fphase = PhaseCalculation(fpeaks);

% jade
W = jadeR(data);
s1 = W*data;

% % % W = jadeR(data2);
% % % s2 = W*data2;
s2 = PiCA(data2,fpeaks);

PlotECG(data,8,'b');
PlotECG(data2,8,'r');

PlotECG(s1,8,'b');
PlotECG(s2,8,'r');

fECGmean1 = zeros(size(data,1),bins);
fECGmean2 = zeros(size(data,1),bins);
fECGsd1 = zeros(size(data,1),bins);
fECGsd2 = zeros(size(data,1),bins);
fmeanphase1 = zeros(size(data,1),bins);
fmeanphase2 = zeros(size(data,1),bins);
for i = 1:size(data,1)
    [fECGmean1(i,:),fECGsd1(i,:),fmeanphase1(i,:)] = MeanECGExtraction(s1(i,:),fphase,bins,1);
    [fECGmean2(i,:),fECGsd2(i,:),fmeanphase2(i,:)] = MeanECGExtraction(s2(i,:),fphase,bins,1);
    
    figure;
    subplot(121);
    errorbar(fECGmean1(i,:),fECGsd1(i,:)/2,'b');
    hold on;
    plot(fECGmean1(i,:),'r','linewidth',2);
    grid on;
    xlabel('Sample index');
    ylabel('Amplitude (mV)');
    axis tight;
    a = axis;
    
    subplot(122);
    errorbar(fECGmean2(i,:),fECGsd2(i,:)/2,'b');
    hold on;
    plot(fECGmean2(i,:),'r','linewidth',2);
    grid on;
    xlabel('Sample index');
    ylabel('Amplitude (mV)');
    %     axis(a);
    axis tight;
end

fECGmean1 = fECGmean1/max(max(abs(fECGmean1)));
fECGmean2 = fECGmean2/max(max(abs(fECGmean2)));

figure;
plot(fECGmean1','linewidth',2);
set(gca,'Box','On','FontSize',16);
xlabel('Sample index','FontSize',16);
ylabel('Normalized Amplitude','FontSize',16);
% % % a = axis;
% % % a(2) = length(fECGmean1);
% % % axis(a);
axis tight;
grid;
set(gcf,'Position',[520 471 560 627]);

figure;
plot(fECGmean2','linewidth',2);
set(gca,'Box','On','FontSize',16);
xlabel('Sample index','FontSize',16);
ylabel('Normalized Amplitude','FontSize',16);
% % % a = axis;
% % % a(2) = length(fECGmean1);
% % % axis(a);
axis tight;
grid;
set(gcf,'Position',[520 471 560 627]);

% % % figure;
% % % hold on;
% % % plot(fetref);
% % % plot(fpeaks.*fetref,'ro');

