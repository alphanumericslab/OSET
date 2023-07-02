function [lambda AIC MDL NEW ENSTh ENS fhand] = EstimateDimension(x, noisevar, flagplot)
% estimates the dimensions of the signal space from noisy data
% Reza Sameni
% May 2011
%
% Reference:
% Raj Rao Nadakuditi and Alan Edelman, "Sample Eigenvalue Based Detection
% of High-Dimensional Signals in White Noise Using Relatively Few Samples"
% IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 56, NO. 7, JULY 2008
%

L = size(x,1);
N = size(x,2);
R = cov(x',1);
[~, d] = eig(R);
lambda = sort(diag(d),'descend');
P = min(N,L);
AIC = zeros(1,P);
MDL = zeros(1,P);
NEW = zeros(1,P);
for k = 0:P-1
    g = geomean(lambda(k+1:L));
    a = mean(lambda(k+1:L));
    AIC(k+1) = -2*(L-k)*N*log(g/a) + 2*k*(2*L - k);
    MDL(k+1) = -(L-k)*N*log(g/a) + 0.5*k*(2*L - k)*log(N);

    beta = 1;
    t = L*((L-k)*sum(lambda(k+1:L).^2)/sum(lambda(k+1:L))^2 - (1 + L/N)) - (2/beta - 1)* L/N;
    NEW(k+1) = beta/4*(N/L)^2*t^2 + 2*(k+1);
end

% effective number of signals threshold:
ENSTh = 1 + sqrt(L/N);

K = find(lambda > noisevar*ENSTh);
ENS = length(K);
fhand = [];
if(flagplot)
    nn = 1:P;
    f1 = figure;
    hold on;
    plot(nn,log(lambda/lambda(1)),'linewidth',3);
    plot(nn,log(lambda/lambda(1)),'ro','linewidth',3);
    xlabel('n');
    ylabel('$\lambda_n$ (dB)','interpreter','latex');
    title('normalized eigenvalues');
    grid

    f2 = figure;
    hold on;
    plot(log(AIC),'-','linewidth',4,'color',0*ones(1,3));
    plot(log(MDL),'-.','linewidth',4,'color',.6*ones(1,3));
    plot(log(NEW),'--','linewidth',4,'color',.3*ones(1,3));
    grid
    legend('AIC','MDL','NE');
    set(gca,'box','on');
    set(gca,'fontsize',14);
    xlabel('m');
    hold off;

    f3 = figure;
    hold on;
    plot(nn,lambda,'linewidth',3,'color',.2*ones(1,3)); % all channels should be normalized to have the same variance
    plot(nn,lambda,'.','markersize',30,'color',.2*ones(1,3)); % all channels should be normalized to have the same variance
    plot(nn,ENSTh(ones(1,length(nn)))+noisevar, '--','linewidth',4,'color',.5*ones(1,3));
    text(round(max(nn)/2),ENSTh + noisevar,'$\delta^2(1 + \sqrt{N_c/T})$','FontSize',18,'interpreter','latex');
    axis('tight');
    grid
    set(gca,'fontsize',16,'box','on');
    set(gca,'yscale','log');
    xlabel('$n$','interpreter','latex');
    ylabel('$\log(l_n)$','interpreter','latex');
    hold off;
    fhand = [f1 f2 f3];
end