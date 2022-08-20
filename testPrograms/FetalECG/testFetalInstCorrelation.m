% Test the instantaneous correlation of different signals
% Reza Sameni, 2009
%

clc;
clear;
close all;


% Typical ECG data. Download from: http://homes.esat.kuleuven.be/~smc/daisy/daisydata.html
load('FOETAL_ECG.dat'); data = FOETAL_ECG(:,2:end)'; time = FOETAL_ECG(:,1)'; clear FOETAL_ECG; fs = 250;

L = size(data,1);
N = size(data,2);

b = LPFilter(data,.7/fs); % Using OSET tools. Download this function and others from: http://ecg.sharif.edu/
dat = data - b;
dat = dat./(std(dat,[],2)*ones(1,N));

% W =  jadeR(dat); % Cardoso's JADE algorithm, from: http://www.tsi.enst.fr/~cardoso/guidesepsou.html
% dat = W*dat;


tau = 100:220;
T = 20;

ch = 1;
r = InstCorr(dat(ch, :),tau);

% r(r < .1) = .1;

R = filter(ones(1,T), T, r, [], 1);
E = sqrt(filter(ones(1,T), T, r.^2, [], 1));
ee = max(E, [], 2);
R = R./(ee*ones(1, size(R, 2)));
% for i = 1:size(r,1),
%     rr = filter(ones(1,T), T, r(i, :), [], 1);
%
% %     R(i, :) = rr./(max()
% end

n = (0:length(dat)-1);

figure
mesh(n/fs,tau/fs,r');
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('t');
ylabel('tau');

figure
mesh(n/fs,tau/fs,R');
% axis xy; axis tight; colormap(jet); view(0,90);
xlabel('t');
ylabel('tau');

% PlotECG(dat,4,'k',fs);

% % % ref = dat(1,:);
% % % f = 1.1;
% % % peaksm = PeakDetection(ref,f/fs);
% % %
% % % W =  jadeR(dat); % Cardoso's JADE algorithm, from: http://www.tsi.enst.fr/~cardoso/guidesepsou.html
% % % ss = W*dat;
% % % ref = ss(8,:);
% % % f = 2;
% % % peaksf = PeakDetection(ref,f/fs);
% % %
% % % [sm,W,A] = PiCA(dat,peaksm);
% % % [sf,W,A] = PiCA(dat,peaksf);
% % %
% % % [T0,T1] = SynchPhaseTimes2(peaksm);
% % %
% % % % % % R0 = dat*dat';
% % % % % % e = exp(1j*2*pi.fs/(T1 - T0);
% % % % % % R1 = dat*(dat.*exp(1j*2*pi*e
% % %
% % %
% % % PlotECG(dat,4,'k',fs);  % Raw data
% % % PlotECG(ss,4,'b',fs);   % ECG decomposed by JADE
% % % PlotECG(sm,4,'r',fs);   % ECG decomposed by PiCA (maternal synchronous)
% % % PlotECG(sf,4,'m',fs);   % ECG decomposed by PiCA (fetal synchronous)
