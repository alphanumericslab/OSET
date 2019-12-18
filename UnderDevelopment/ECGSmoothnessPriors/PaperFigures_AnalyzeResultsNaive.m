close all
clc
clear
% results = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB04_NaiveSweepGammaBased_PartialResults_Backup.txt');
% results = [zeros(size(results.data(:,1),1), 2) results.data(:,1:end)];
% results = load('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia04_NaiveSweepGammaBased_Backup.txt');
results = load('J:\ECGData\Normal\AllDataSmootherResultsNormal04_NaiveSweepGammaBased_Backup.txt');

snr_list = 30 : -3: -15;

snr = results(:, 8);
order = results(:, 5);
snr0 = results(:, 9);
snrimp1 = results(:, 10);
optparam = results(:, 11);

snrimp1_mean = zeros(1, length(snr_list));
snrimp1_std = zeros(1, length(snr_list));
optparam1_mean = zeros(1, length(snr_list));
optparam1_std = zeros(1, length(snr_list));
snrimp2_mean = zeros(1, length(snr_list));
snrimp2_std = zeros(1, length(snr_list));
optparam2_mean = zeros(1, length(snr_list));
optparam2_std = zeros(1, length(snr_list));
snrimp3_mean = zeros(1, length(snr_list));
snrimp3_std = zeros(1, length(snr_list));
optparam3_mean = zeros(1, length(snr_list));
optparam3_std = zeros(1, length(snr_list));
for i = 1 : length(snr_list),
    I1 = find(snr == snr_list(i) & order == 2);
    snrimp1_mean(i) = mean(snrimp1(I1));
    snrimp1_std(i) = std(snrimp1(I1));
    optparam1_mean(i) = median(optparam(I1));
    optparam1_std(i) = std(optparam(I1));

    I2 = find(snr == snr_list(i) & order == 4);
    snrimp2_mean(i) = mean(snrimp1(I2));
    snrimp2_std(i) = std(snrimp1(I2));
    optparam2_mean(i) = median(optparam(I2));
    optparam2_std(i) = std(optparam(I2));

    I3 = find(snr == snr_list(i) & order == 6);
    snrimp3_mean(i) = mean(snrimp1(I3));
    snrimp3_std(i) = std(snrimp1(I3));
    optparam3_mean(i) = median(optparam(I3));
    optparam3_std(i) = std(optparam(I3));
end

figure
hold on
% errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
plot(snr_list, snrimp1_mean, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
plot(snr_list, snrimp2_mean, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
plot(snr_list, snrimp3_mean, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
grid
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
axis tight
a = axis;
a(1) = -16;%snr_list(end);
a(2) = 31;%snr_list(1);
% a(3) = 0;
% a(4) = 8;
axis(a);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR Mean (dB)');

figure
hold on
plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
plot(snr_list, snrimp2_std, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
plot(snr_list, snrimp3_std, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
grid
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
set(gcf, 'Position', [520 378 560 219])
axis tight
a = axis;
a(1) = -16;%snr_list(end);
a(2) = 31;%snr_list(1);
% a(3) = 0;
% a(4) = 8;
axis(a);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR SD (dB)');

figure
hold on
% errorbar(snr_list, optparam1_mean, optparam1_std, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% errorbar(snr_list, optparam2_mean, optparam2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% errorbar(snr_list, optparam3_mean, optparam3_std, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
plot(snr_list, optparam1_mean, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
plot(snr_list, optparam2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
plot(snr_list, optparam3_mean, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
set(gca, 'yscale', 'log');
grid
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
set(gcf, 'Position', [520 81 560 219])
axis tight
a = axis;
a(1) = -16;%snr_list(end);
a(2) = 31;%snr_list(1);
% a(3) = 0;
% a(4) = 8;
axis(a);
xlabel('Input SNR (dB)');
ylabel('\gamma_0^*');

