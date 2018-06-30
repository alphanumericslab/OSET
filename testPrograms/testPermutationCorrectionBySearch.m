% test algorithm for ICA permutation compensation by search
% UNDER TEST
% Reza Sameni, Copyright 2015
%

clear;
close all;
clc;

N = 1000;
L = 28;
itr = 1;
x1 = randn(L, N);
x2 = randn(L, N);

% load FetalECGDaISy
load('FOETAL_ECG.dat'); data = FOETAL_ECG(:, 2:end)'; fs = 250;
x10 = fastica(data(:, 1:end/2),'approach', 'symm', 'displayMode', 'off');
x20 = fastica(data(:, end/2+1:end),'approach', 'symm', 'displayMode', 'off');
x1 = x10;
x2 = x20;
N = size(x1, 2);
L = size(x1, 1);

Cx0 = x1*x2'/N;
% cost_opt = (sum(Cx0(:)) - trace(Cx0))/sum(Cx0(:));
cost_opt = trace(Cx0);
I_opt = 1:L;
P = eye(L);
Cx = Cx0;
for k = 1:itr,
    I = randperm(L);
    xtmp = P(:, I)*x2;
    Cxtmp = x1*xtmp'/N;
    %     ccost = (sum(Cx(:)) - trace(Cx))/sum(Cx(:));
    ccost = trace(Cxtmp);
    if( ccost > cost_opt) % && ccost > 0
        x2 = xtmp;
        Cx = Cxtmp;
        %         I_opt = I;
        cost_opt = ccost
        
        CCx0 = 32*(Cx0 - min(Cx0(:)))/max(Cx0(:));
        CCx = 32*(Cx - min(Cx(:)))/max(Cx(:));
        figure
        subplot(121);
        image(CCx0');
        axis square
        colormap gray
        subplot(122);
        image(CCx');
        axis square
        colormap gray
    end
end
% Cy = P(I_opt, :)*Cx;

[U S V] = svd(Cx0);
Q = V*U';
Cy = Q*Cx0;

Cz = Cx0;
for k = 1 : 10,
    for i = 1:size(Cx0, 1)
        temp = Cz;
        [v, I] = max(temp(i, :));
        temp(:, [i I]) = temp(:, [I i]);
        if(trace(temp) > trace(Cz))
            Cz(:, [i I]) = Cz(:, [I i]);
        end
    end
end

CCx0 = 32*(Cx0 - min(Cx0(:)))/max(Cx0(:));
CCx = 32*(Cx - min(Cx(:)))/max(Cx(:));
CCy = 32*(Cy - min(Cy(:)))/max(Cy(:));
CCz = 32*(Cz - min(Cz(:)))/max(Cz(:));

figure
subplot(121);
image(CCx0');
axis square
colormap gray
subplot(122);
image(CCx');
axis square
colormap gray

figure
subplot(121);
image(CCx0');
axis square
colormap gray
subplot(122);
image(CCy');
axis square
colormap gray

figure
subplot(121);
image(CCx0');
axis square
colormap gray
subplot(122);
image(CCz');
axis square
colormap gray

% PlotECG(x10, 8, 'b', fs);
% PlotECG(x20, 8, 'r', fs);
%
% PlotECG(x10, 8, 'b', fs);
% PlotECG(Q*x20, 8, 'm', fs);