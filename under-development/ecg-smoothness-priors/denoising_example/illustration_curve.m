% A script for generating illustrations in the following publication: 
%   Sameni, R. (2017). Online filtering using piecewise smoothness priors: 
%       Application to normal and abnormal electrocardiogram denoising. In 
%       Signal Processing (Vol. 133, pp. 52–63). Elsevier BV.
%       https://doi.org/10.1016/j.sigpro.2016.10.019


close all
clc
clear
randn('seed', 0);

t = 0.2:.003:1.5;
s = (t-0.5).^3 - 1*(t-0.55).^2;

x = s + .015*randn(size(s));

figure
hold on
plot(t, x, 'color', 'r');%0.5*ones(1,3));
plot(t, s, 'linewidth', 3, 'color', 0.25*ones(1,3));
axis 'tight'
set(gca, 'fontsize', 14);
h = line(0.5*[1 1], [-1 1]);set(h, 'linestyle', '--', 'linewidth', 3, 'color', 'b');%0.5*ones(1,3));
h = line(1.2*[1 1], [-1 1]);set(h, 'linestyle', '--', 'linewidth', 3, 'color', 'b');%0.5*ones(1,3));
h = line(0.815*[1 1], [-1 1]);set(h, 'linestyle', '--', 'linewidth', 3, 'color', 0.8*ones(1,3));
h = line(0.35*[1 1], [-1 1]);set(h, 'linestyle', '--', 'linewidth', 2, 'color', 0.7*ones(1,3));
h = line(1.35*[1 1], [-1 1]);set(h, 'linestyle', '--', 'linewidth', 2, 'color', 0.7*ones(1,3));

% replot the figures over the lines
plot(t, x, 'color', 'r');%0.5*ones(1,3));
plot(t, s, 'linewidth', 3, 'color', 0.25*ones(1,3));

I1 = find(t >= .35, 1);
I2 = find(t >= .5, 1);
I3 = find(t >= 1.2, 1);
I4 = find(t >= 1.35, 1);
plot(t(I1:I2), s(I1:I2), 'linewidth', 3, 'color', 0.6*ones(1,3));
plot(t(I3:I4), s(I3:I4), 'linewidth', 3, 'color', 0.6*ones(1,3));

% text(0.3, 0.091, 'S_{k-1}', 'fontsize', 14);
% text(0.82, 0.091, 'S_{k}', 'fontsize', 14);
% text(1.3, 0.091, 'S_{k+1}', 'fontsize', 14);
set(gca, 'XTickMode', 'manual');
set(gca, 'XTickLabelMode', 'manual');
set(gca, 'XTick', [])
set(gca, 'XTickLabel', {''});
set(gca, 'YTickMode', 'manual');
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTick', [])
set(gca, 'YTickLabel', {''});
% set(gca, 'XTick', [.4, .5, 1.2, 1.3])
% set(gca, 'XTickLabel', {'|', 'pk', 'p_{k+1}', '|', '|'});
% text(.3, -.175, 'S_{k-1}', 'fontsize', 14);
% text(.82, -.175, 'S_{k}', 'fontsize', 14);
% text(1.32, -.175, 'S_{k+1}', 'fontsize', 14);

text(.4, .073, '$q$', 'fontsize', 14, 'interpreter', 'latex');
[Xf1, Yf1] = ds2nfu(.35, .065);
[Xf2, Yf2] = ds2nfu(.5, .065);
annotation('line',[Xf1, Xf2],[Yf1, Yf2]);

text(1.25, .073, '$q$', 'fontsize', 14, 'interpreter', 'latex');
[Xf1, Yf1] = ds2nfu(1.2, .065);
[Xf2, Yf2] = ds2nfu(1.35, .065);
annotation('line',[Xf1, Xf2],[Yf1, Yf2]);

text(.25, .02, '$\mathbf{\alpha}_{k}$', 'fontsize', 14, 'interpreter', 'latex');
[Xf1, Yf1] = ds2nfu(.28, .015);
[Xf2, Yf2] = ds2nfu(.4, -.02);
annotation('arrow',[Xf1, Xf2],[Yf1, Yf2]);

text(1.4, -.09, '$\mathbf{\beta}_{k}$', 'fontsize', 14, 'interpreter', 'latex');
[Xf1, Yf1] = ds2nfu(1.38, -.09);
[Xf2, Yf2] = ds2nfu(1.3, -.05);
annotation('arrow',[Xf1, Xf2],[Yf1, Yf2]);

text(.92, .01, '$\mathbf{x}_{k}$', 'fontsize', 14, 'interpreter', 'latex');
[Xf1, Yf1] = ds2nfu(.95, 0);
[Xf2, Yf2] = ds2nfu(.97, -.04);
annotation('arrow',[Xf1, Xf2],[Yf1, Yf2]);

text(.6, -.07, '$\mathbf{\theta}_{k}$', 'fontsize', 14, 'interpreter', 'latex');
[Xf1, Yf1] = ds2nfu(.63, -.058);
[Xf2, Yf2] = ds2nfu(.7, -.01);
annotation('arrow',[Xf1, Xf2],[Yf1, Yf2]);

plot([t(I1), t(I2), t(I3), t(I4)], [s(I1), s(I2), s(I3), s(I4)], 'ko', 'MarkerSize', 8, 'markerfacecolor', 'k')

text(.2, -.175, '$S_{k-1}$', 'fontsize', 14, 'interpreter', 'latex');
text(.84, -.175, '$S_{k}$', 'fontsize', 14, 'interpreter', 'latex');
text(1.4, -.175, '$S_{k+1}$', 'fontsize', 14, 'interpreter', 'latex');
text(.48, -.175, '$p_{k}$', 'fontsize', 14, 'interpreter', 'latex');
% text(.32, -.175, '$p_{k}-q$', 'fontsize', 14, 'interpreter', 'latex');
text(1.18, -.175, '$p_{k+1}$', 'fontsize', 14, 'interpreter', 'latex');

text(.48, .13, '$\tilde{S}_{k}$', 'fontsize', 14, 'interpreter', 'latex');
text(1.15, .13, '$\tilde{S}_{k+1}$', 'fontsize', 14, 'interpreter', 'latex');
text(0.78, .13, '$\tilde{p}_{k+1}$', 'fontsize', 14, 'interpreter', 'latex');

[Xf1, Yf1] = ds2nfu(.205, .112);
[Xf2, Yf2] = ds2nfu(.81, .112);
annotation('doublearrow',[Xf1, Xf2],[Yf1, Yf2]);

[Xf1, Yf1] = ds2nfu(.82, .112);
[Xf2, Yf2] = ds2nfu(1.5, .112);
annotation('doublearrow',[Xf1, Xf2],[Yf1, Yf2]);

[Xf1, Yf1] = ds2nfu(.5, -.15);
[Xf2, Yf2] = ds2nfu(1.2, -.15);
annotation('doublearrow',[Xf1, Xf2],[Yf1, Yf2]);

[Xf1, Yf1] = ds2nfu(1.2, -.15);
[Xf2, Yf2] = ds2nfu(1.5, -.15);
annotation('arrow',[Xf2, Xf1],[Yf2, Yf1]);

[Xf1, Yf1] = ds2nfu(0.2, -.15);
[Xf2, Yf2] = ds2nfu(.5, -.15);
annotation('arrow',[Xf1, Xf2],[Yf1, Yf2]);

axis([0.2000    1.4990   -0.16    0.1145])
