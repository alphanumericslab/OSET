% test algorithm for ICA permutation compensation by SVD
% UNDER TEST
% Reza Sameni, Copyright 2015
%

clear;
close all;
clc;

% N = 1000;
% L = 8;
% itr = 1;
% x1 = randn(L, N);
% x2 = randn(L, N);
% fs = 250;

% load('FetalECGDaISy.mat');
load('FOETAL_ECG.dat'); data = FOETAL_ECG(:, 2:end)'; fs = 250;

x1 = fastica(data(:, 1:end/2),'approach', 'symm', 'displayMode', 'off');
x2 = fastica(data(:, end/2+1:end),'approach', 'symm', 'displayMode', 'off');
N = size(x1, 2);
L = size(x1, 1);

Cx1 = x1*x1'/N;
Cx2 = x2*x2'/N;

[V1, D1] = eig(Cx1);
[V2, D2] = eig(Cx2);

s1 = x1;%diag(diag(D1).^-0.5)*V1'*x1;
s2 = x2;%diag(diag(D2).^-0.5)*V2'*x2;

% s1*s1'/N
% s2*s2'/N

C12 = s1*s2'/N;
[U S V] = svd(C12);
Q = V*U';
s22 = Q'*s2;

C1 = s1*s2'/N;
C2 = s1*s22'/N;

PlotECG(s1, 8, 'b', fs);
PlotECG(s2, 8, 'r', fs);
PlotECG(s22, 8, 'm', fs);

CC1 = 32*(C1 - min(C1(:)))/max(C1(:));
CC2 = 32*(C2 - min(C2(:)))/max(C2(:));

figure
subplot(121);
% image(CC1');
mesh(C1');
axis square
colormap gray
subplot(122);
image(CC2');
% mesh(C2');
axis square
colormap gray

% % cost_opt = (sum(Cx0(:)) - trace(Cx0))/sum(Cx0(:));
% cost_opt = trace(Cx0);
% I_opt = 1:L;
% P = eye(L);
% Cx = Cx0;
% for k = 1:itr,
%     I = randperm(L);
%     xtmp = P(:, I)*x2;
%     Cxtmp = x1*xtmp'/N;
%     %     ccost = (sum(Cx(:)) - trace(Cx))/sum(Cx(:));
%     ccost = trace(Cxtmp);
%     if( ccost > cost_opt) % && ccost > 0
%         x2 = xtmp;
%         Cx = Cxtmp;
%         %         I_opt = I;
%         cost_opt = ccost
%         
%         CCx0 = 32*(Cx0 - min(Cx0(:)))/max(Cx0(:));
%         CCx = 32*(Cx - min(Cx(:)))/max(Cx(:));
%         figure
%         subplot(121);
%         image(CCx0');
%         axis square
%         colormap gray
%         subplot(122);
%         image(CCx');
%         axis square
%         colormap gray
%     end
% end
% % Cy = P(I_opt, :)*Cx;
% 
% [U S V] = svd(Cx0);
% Q = V*U';
% Cy = Q*Cx0;
% 
% Cz = Cx0;
% for k = 1 : 10,
%     for i = 1:size(Cx0, 1)
%         temp = Cz;
%         [~, I] = max(temp(i, :));
%         temp(:, [i I]) = temp(:, [I i]);
%         if(trace(temp) > trace(Cz))
%             Cz(:, [i I]) = Cz(:, [I i]);
%         end
%     end
% end
% 
% CCx0 = 32*(Cx0 - min(Cx0(:)))/max(Cx0(:));
% CCx = 32*(Cx - min(Cx(:)))/max(Cx(:));
% CCy = 32*(Cy - min(Cy(:)))/max(Cy(:));
% CCz = 32*(Cz - min(Cz(:)))/max(Cz(:));
% 
% figure
% subplot(121);
% image(CCx0');
% axis square
% colormap gray
% subplot(122);
% image(CCx');
% axis square
% colormap gray
% 
% figure
% subplot(121);
% image(CCx0');
% axis square
% colormap gray
% subplot(122);
% image(CCy');
% axis square
% colormap gray
% 
% figure
% subplot(121);
% image(CCx0');
% axis square
% colormap gray
% subplot(122);
% image(CCz');
% axis square
% colormap gray
% 
% % PlotECG(x10, 8, 'b', fs);
% % PlotECG(x20, 8, 'r', fs);
% %
% % PlotECG(x10, 8, 'b', fs);
% PlotECG(Q*x20, 8, 'm', fs);