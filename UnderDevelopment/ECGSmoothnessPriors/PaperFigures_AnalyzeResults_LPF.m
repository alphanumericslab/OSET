close all
clc
clear
results = load('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia05_ComparisonWithLPF_Backup.txt');
% results = load('J:\ECGData\Normal\AllDataSmootherResultsNormal05_ComparisonWithLPF_Backup.txt');
% results = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB05_ComparisonWithLPF_Backup.txt');
% results = [zeros(size(results.data(:,1),1), 2) results.data(:,1:end)];

% snr_list = 30 : -3: -12;
snr_list = 30 : -3: -15;

snr = results(:, 8);
snr0 = results(:, 9);
snrimp1 = results(:, 10);
% snrimp2 = results(:, 10);
% snrimp3 = results(:, 11);
optparam = results(:, 11);

snrimp1_mean = zeros(1, length(snr_list));
% snrimp2_mean = zeros(1, length(snr_list));
% snrimp3_mean = zeros(1, length(snr_list));
snrimp1_std = zeros(1, length(snr_list));
% snrimp2_std = zeros(1, length(snr_list));
% snrimp3_std = zeros(1, length(snr_list));
optparam_mean = zeros(1, length(snr_list));
optparam_std = zeros(1, length(snr_list));
for i = 1 : length(snr_list),
    I = find(snr == snr_list(i));
    snrimp1_mean(i) = mean(snrimp1(I));
%     snrimp2_mean(i) = mean(snrimp2(I));
%     snrimp3_mean(i) = mean(snrimp3(I));
    snrimp1_std(i) = std(snrimp1(I));
%     snrimp2_std(i) = std(snrimp2(I));
%     snrimp3_std(i) = std(snrimp3(I));
    optparam_mean(i) = mean(optparam(I));
    optparam_std(i) = std(optparam(I));
end

figure
hold on
% errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
plot(snr_list, snrimp1_mean, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
grid
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
a = axis;
a(1) = -15;%snr_list(end);
a(2) = 33;%snr_list(1);
a(3) = 0;
a(4) = 8;
% axis(a);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR Mean (dB)');

figure
hold on
plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
grid
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
a = axis;
a(1) = -15;%snr_list(end);
a(2) = 33;%snr_list(1);
a(3) = 0;
a(4) = 2;
% axis(a);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR SD (dB)');

figure
hold on
errorbar(snr_list, optparam_mean, optparam_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
grid
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
% a = axis;
% a(1) = -15;%snr_list(end);
% a(2) = 33;%snr_list(1);
% a(3) = 0;
% a(4) = 2;
% axis(a);
xlabel('Input SNR (dB)');
ylabel('Optim Param');
