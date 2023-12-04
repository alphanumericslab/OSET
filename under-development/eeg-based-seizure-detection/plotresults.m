% clear
close all;
% r = load('results_testProcessSeizureEEG10_Multiple_ICA_all_cases.txt');
% r = load('results_testProcessSeizureEEG10_Multiple_ICA_all_cases_with_angles.txt');
% r = load('results_ICA.txt');

r = data;
% bins = 0:.01:1;
% bins = -10:.5:50;
bins = 30;%0:1:100;
column = 4;

subject = 1;%(r(:, 1) == 5);
mode1 = (r(:, 3) == 1); % interictal
mode2 = (r(:, 3) == 2); % preictal
mode3 = (r(:, 3) == 3); % test

[n1, r1] = hist(r(mode1 & subject, column), bins);
[n2, r2] = hist(r(mode2 & subject, column), bins);
[n3, r3] = hist(r(mode3 & subject, column), bins);

N1 = sum(n1);
N2 = sum(n2);
N3 = sum(n3);

figure
hold on
bar(r1, n1/N1);
bar(r2, n2/N2, 'r');
% bar(r3, n3/N3, 'g');

% h = findobj(gca,'Type','patch');
% set(h(1), 'FaceColor', 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'EdgeAlpha', 0.9);
% set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.9);
% set(h(3), 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.9);
% legend('interictal', 'preictal', 'test');
grid
