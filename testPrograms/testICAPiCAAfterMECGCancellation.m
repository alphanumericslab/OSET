% applying ICA and Pi-CA after maternal ECG removal

clear;
close all;

% load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; time = FOETAL_ECG(:,1)'; clear FOETAL_ECG; fs = 250;
load('data.txt');
time = data(:,1)';
data = data(:,2:end)';
fs = 1000;
data = data + 32768;
data(data>32768) = data(data>32768) - 65536;

L1 = size(data,1);
L2 = size(data,2);
bins = fs;

% R-peak detection
% ref = data(5,:);
ref = data(9,:);
b = LPFilter(ref,1.5/fs);
ref = ref - b;
f = 1.7;
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

for i = 1:Itr
    
    % periodic component analysis stage
    A = dat(:,T0)*dat(:,T1)';
    B = dat(:,T0)*dat(:,T0)';
    
    A = (A+A')/2;
    B = (B+B')/2;
    
    STOPCRITERION = trace(A)/trace(B);
    
    [V,D] = eig(A,B,'chol');
    
    d = diag(D);
    [Y,I] = sort(d);
    I = I(end:-1:1);
    
    W = V(:,I)';
    dat = W*dat;
    %     PlotECG(dat,8,'b');
    
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
    
    % % %     ekfsnrimp(itr) = 10*log10(NoisePower/mean((noise-resEKF).^2)) - CCsnr(j);
    % % %     ekssnrimp(itr) = 10*log10(NoisePower/mean((noise-resEKS).^2)) - CCsnr(j);
    % % %
    % % %     crcoef = corrcoef(xxx(time1),xxx(time2));                                       inputPM(itr) = crcoef(1,2);
    % % %     crcoef = corrcoef(resEKF(time1)-bsline(time1)',resEKF(time2)-bsline(time2)');   ekfoutputPM(itr) = crcoef(1,2);
    % % %     crcoef = corrcoef(resEKS(time1)-bsline(time1)',resEKS(time2)-bsline(time2)');   eksoutputPM(itr) = crcoef(1,2);
    
end

data2 = wholedata(:,:,end);

fetref = data2(1,:);
fpeaks = PeakDetection(fetref,2.2/fs);
fphase = PhaseCalculation(fpeaks);

% jade
W = jadeR(data2);
s1 = W*data2;

% % % W = jadeR(data2);
% % % s2 = W*data2;
s2 = PiCA(data2,fpeaks);

PlotECG(data,4,'b',fs);
PlotECG(data2,4,'r',fs);

PlotECG(s1,4,'b',fs);
PlotECG(s2,4,'m',fs);

fECGmean1 = zeros(size(data,1),bins);
fECGmean2 = zeros(size(data,1),bins);
fECGsd1 = zeros(size(data,1),bins);
fECGsd2 = zeros(size(data,1),bins);
fmeanphase1 = zeros(size(data,1),bins);
fmeanphase2 = zeros(size(data,1),bins);


% % % figure;
% % % for i = 1:L1,
% % %     subplot(L1,1,i);
% % %     plot(time,s1,'k','linewidth',1);
% % %     grid;
% % %     set(gca,'Box','On','FontSize',16);
% % %     if (i<L1)
% % %         set(gca,'XTickLabel',[]);
% % %     end
% % %     %     ylabel(['IC_',num2str(i)],'FontSize',16);
% % %     if (itr ==1 && i==1),
% % %         title('Original Data');
% % %     end
% % %     if (itr ==Itr+1 && i==1),
% % %         title('After Maternal Artifact Removal');
% % %     end
% % % end
% % % xlabel('Time (s)','FontSize',16);
% % %
% % % figure;
% % % for i = 1:L1,
% % %     subplot(L1,1,i);
% % %     plot(time,s2,'k','linewidth',1);
% % %     grid;
% % %     set(gca,'Box','On','FontSize',16);
% % %     if (i<L1)
% % %         set(gca,'XTickLabel',[]);
% % %     end
% % %     %     ylabel(['IC_',num2str(i)],'FontSize',16);
% % %     if (itr ==1 && i==1),
% % %         title('Original Data');
% % %     end
% % %     if (itr ==Itr+1 && i==1),
% % %         title('After Maternal Artifact Removal');
% % %     end
% % % end
% % % xlabel('Time (s)','FontSize',16);
