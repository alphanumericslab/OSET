% test the Spikes and Bumps phenomenon in ICA
clear;
close all;
randn('seed',0);
rand('seed',0);

x = cumsum(randn(40,1000),2);

[s, A, W] = fastica(x,'approach', 'symm','displayMode', 'off','epsilon',.001);
% % % W =  jadeR(d); s = real(W*d); A = pinv(W);

PlotECG(s',4,'b');

figure;
plot(s(14,:),'lineWidth',2);
grid;
xlabel('Samples','FontSize',16);
ylabel('Amplitude','FontSize',16);
set(gca,'Box','On','FontSize',16);
set(gcf,'Position',[438 350 1200 420]);
