close all
clc
clear
% results1 = load('J:\ECGData\PTB\AllDataSmootherResultsPTB07_NaiveSweepGammaBasedOptimized_Backup.txt');
% results1 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal04_NaiveSweepGammaBased_Backup.txt');
% results1 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal07_NaiveSweepGammaBasedOptimized.txt');

% results = load('J:\ECGData\Normal\AllDataSmootherResultsNormal08_KnownVarianceOptimized_Backup.txt');
results = load('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia08_KnownVarianceOptimized_Backup.txt');

snr_list = 30 : -3: -15;
diff_h = [2 4 6];
wlen = [100e-3 200e-3];%%%100e-3; % window length (s)
trial_num = 1;
NoisePowerOverEstimationFactor = [0.8 1 1.2];%%%sqrt(fs*wlen)];
ForgettingFactor = 1;%%%0.9 : 0.1 : 1; % forgetting rate of the smoothness factor calculated in each block

wlens = results(:, 7);
orders = results(:, 8);
forgettingfactors = results(:, 9);
overestim = results(:, 10);
snr = results(:, 11);
snr0 = results(:, 12);
snrimp1 = results(:, 13);
snrimp2 = results(:, 14);
optparam1 = results(:, 15);
optparam2 = results(:, 16);

snrimp1_mean = zeros(length(snr_list), length(wlen), length(diff_h), length(ForgettingFactor), length(NoisePowerOverEstimationFactor));
snrimp1_std = zeros(length(snr_list), length(wlen), length(diff_h), length(ForgettingFactor), length(NoisePowerOverEstimationFactor));
snrimp2_mean = zeros(length(snr_list), length(wlen), length(diff_h), length(ForgettingFactor), length(NoisePowerOverEstimationFactor));
snrimp2_std = zeros(length(snr_list), length(wlen), length(diff_h), length(ForgettingFactor), length(NoisePowerOverEstimationFactor));
% snrimp_sum = zeros(length(snr_list), length(wlens), length(diff_h), length(ForgettingFactor), length(NoisePowerOverEstimationFactor));

% % % figure
% % % hold on;
% for i = 1 : length(snrimp1)
for ss = 1 : length(snr_list),
    for ww = 1 : length(wlen),
        for oo = 1 : length(diff_h),
            for ff = 1 : length(ForgettingFactor),
                for ovr = 1 : length(NoisePowerOverEstimationFactor),
                    I = find(snr == snr_list(ss) & orders == diff_h(oo) & wlens == ww & overestim == ovr & forgettingfactors == ff);
                    snrimp1_mean(ss, ww, oo, ff, ovr) = mean(snrimp1(I));
                    snrimp2_mean(ss, ww, oo, ff, ovr) = mean(snrimp2(I));
                    % % %                     snrimp1_std(ss, ww, oo, ff, ovr) = nanstd(snrimp1(I));
                    %                     optparam1_mean(ss, ww, oo, ff, ovr) = median(optparam1(I));
                    %                     optparam1_std(ss, ww, oo, ff, ovr) = std(optparam1(I));
                    
                    % % %                     snrimp2_mean(ss, ww, oo, ff, ovr) = nanmean(snrimp2(I));
                    % % %                     snrimp2_std(ss, ww, oo, ff, ovr) = nanstd(snrimp2(I));
                    %                     optparam2_mean(ss, ww, oo, ff, ovr) = median(optparam2(I));
                    %                     optparam2_std(ss, ww, oo, ff, ovr) = std(optparam2(I));
                    %                     plot(snr_list, squeeze(snrimp1_mean(:, ww, oo, ff, ovr)), 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
                    %                     plot(snr_list, squeeze(snrimp2_mean(:, ww, oo, ff, ovr)), 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
                end
            end
        end
    end
end

% end
figure
hold on
for ww = 1 : length(wlen),
    for oo = 1 : length(diff_h),
        for ff = 1 : length(ForgettingFactor),
            for ovr = 1 : length(NoisePowerOverEstimationFactor),
                plot(snr_list, squeeze(snrimp1_mean(:, ww, oo, ff, ovr)), 'b', 'marker', 'square', 'LineWidth', 1, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
                plot(snr_list, squeeze(snrimp2_mean(:, ww, oo, ff, ovr)), 'r', 'marker', 'square', 'LineWidth', 1, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
            end
        end
    end
end
grid

snrimp_mean_averaged = mean(snrimp2_mean, 1);
[YY II] = sort(snrimp_mean_averaged(:), 'descend');
[ww, oo, ff, ovr] = ind2sub(size(snrimp_mean_averaged), II(1));
plot(snr_list, squeeze(snrimp2_mean(:, ww, oo, ff, ovr)), 'k', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);

wlen(ww)
diff_h(oo)
ForgettingFactor(ff)
NoisePowerOverEstimationFactor(ovr)

% % % plot(snr_list, squeeze(snrimp1_mean(:, ww, oo, ff, ovr)), 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % plot(snr_list, squeeze(snrimp2_mean(:, ww, oo, ff, ovr)), 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
% % %                     plot(snr_list, snrimp3_mean, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);


% % %                     snrimp1_mean = zeros(1, length(snr_list));
% % %                     snrimp1_std = zeros(1, length(snr_list));
% % %                     optparam1_mean = zeros(1, length(snr_list));
% % %                     optparam1_std = zeros(1, length(snr_list));
% % %                     snrimp2_mean = zeros(1, length(snr_list));
% % %                     snrimp2_std = zeros(1, length(snr_list));
% % %                     optparam2_mean = zeros(1, length(snr_list));
% % %                     optparam2_std = zeros(1, length(snr_list));
% % %                     snrimp3_mean = zeros(1, length(snr_list));
% % %                     snrimp3_std = zeros(1, length(snr_list));
% % %                     optparam3_mean = zeros(1, length(snr_list));
% % %                     optparam3_std = zeros(1, length(snr_list));
% % %                     snrimp4_mean = zeros(1, length(snr_list));
% % %                     snrimp4_std = zeros(1, length(snr_list));
% % %                     optparam4_mean = zeros(1, length(snr_list));
% % %                     optparam4_std = zeros(1, length(snr_list));
% % %                     snrimp5_mean = zeros(1, length(snr_list));
% % %                     snrimp5_std = zeros(1, length(snr_list));
% % %                     optparam5_mean = zeros(1, length(snr_list));
% % %                     optparam5_std = zeros(1, length(snr_list));
% % %                     snrimp6_mean = zeros(1, length(snr_list));
% % %                     snrimp6_std = zeros(1, length(snr_list));
% % %                     optparam6_mean = zeros(1, length(snr_list));
% % %                     optparam6_std = zeros(1, length(snr_list));
% % %                     for i = 1 : length(snr_list),
% % %                         I1 = find(snr == snr_list(i) & order == 2 & wlen == 1 & ovr==2);
% % %                         snrimp1_mean(i) = mean(snrimp2(I1));
% % %                         snrimp1_std(i) = std(snrimp2(I1));
% % %                         optparam1_mean(i) = median(optparam2(I1));
% % %                         optparam1_std(i) = std(optparam2(I1));
% % %
% % %                         I2 = find(snr == snr_list(i) & order == 4 & wlen == 1 & ovr==2);
% % %                         snrimp2_mean(i) = mean(snrimp2(I2));
% % %                         snrimp2_std(i) = std(snrimp2(I2));
% % %                         optparam2_mean(i) = median(optparam2(I2));
% % %                         optparam2_std(i) = std(optparam2(I2));
% % %
% % %                         I3 = find(snr == snr_list(i) & order == 6 & wlen == 1 & ovr==2);
% % %                         snrimp3_mean(i) = mean(snrimp2(I3));
% % %                         snrimp3_std(i) = std(snrimp2(I3));
% % %                         optparam3_mean(i) = median(optparam2(I3));
% % %                         optparam3_std(i) = std(optparam2(I3));
% % %
% % %                         I4 = find(snr == snr_list(i) & order == 2 & wlen == 1 & ovr==1);
% % %                         snrimp4_mean(i) = mean(snrimp2(I4));
% % %                         snrimp4_std(i) = std(snrimp2(I4));
% % %                         optparam4_mean(i) = median(optparam2(I4));
% % %                         optparam4_std(i) = std(optparam2(I4));
% % %
% % %                         I5 = find(snr == snr_list(i) & order == 4 & wlen == 1 & ovr==1);
% % %                         snrimp5_mean(i) = mean(snrimp2(I5));
% % %                         snrimp5_std(i) = std(snrimp2(I5));
% % %                         optparam5_mean(i) = median(optparam2(I5));
% % %                         optparam5_std(i) = std(optparam2(I5));
% % %
% % %                         I6 = find(snr == snr_list(i) & order == 6 & wlen == 1 & ovr==1);
% % %                         snrimp6_mean(i) = mean(snrimp2(I6));
% % %                         snrimp6_std(i) = std(snrimp2(I6));
% % %                         optparam6_mean(i) = median(optparam2(I6));
% % %                         optparam6_std(i) = std(optparam2(I6));
% % %                     end
% % %
% % %                     figure
% % %                     hold on
% % %                     % errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% % %                     % errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% % %                     % plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% % %                     % plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
% % %                     plot(snr_list, snrimp1_mean, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp2_mean, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp3_mean, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp4_mean, 'b--', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp5_mean, 'r--', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp6_mean, 'g--', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     % plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % %                     grid
% % %                     set(gca, 'fontsize', 16);
% % %                     set(gca, 'box', 'on');
% % %                     axis tight
% % %                     a = axis;
% % %                     a(1) = -16;%snr_list(end);
% % %                     a(2) = 31;%snr_list(1);
% % %                     % a(3) = 0;
% % %                     % a(4) = 8;
% % %                     axis(a);
% % %                     xlabel('Input SNR (dB)');
% % %                     ylabel('\Delta SNR Mean (dB)');
% % %
% % %                     figure
% % %                     hold on
% % %                     plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp2_std, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     plot(snr_list, snrimp3_std, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % %                     % plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % %                     grid
% % %                     set(gca, 'fontsize', 16);
% % %                     set(gca, 'box', 'on');
% % %                     set(gcf, 'Position', [520 378 560 219])
% % %                     axis tight
% % %                     a = axis;
% % %                     a(1) = -16;%snr_list(end);
% % %                     a(2) = 31;%snr_list(1);
% % %                     % a(3) = 0;
% % %                     % a(4) = 8;
% % %                     axis(a);
% % %                     xlabel('Input SNR (dB)');
% % %                     ylabel('\Delta SNR SD (dB)');
% % %
% % %                     figure
% % %                     hold on
% % %                     % errorbar(snr_list, optparam1_mean, optparam1_std, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % %                     % errorbar(snr_list, optparam2_mean, optparam2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % %                     % errorbar(snr_list, optparam3_mean, optparam3_std, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % %                     plot(snr_list, optparam1_mean, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % %                     plot(snr_list, optparam2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % %                     plot(snr_list, optparam3_mean, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % %                     set(gca, 'yscale', 'log');
% % %                     grid
% % %                     set(gca, 'fontsize', 16);
% % %                     set(gca, 'box', 'on');
% % %                     set(gcf, 'Position', [520 81 560 219])
% % %                     axis tight
% % %                     a = axis;
% % %                     a(1) = -16;%snr_list(end);
% % %                     a(2) = 31;%snr_list(1);
% % %                     % a(3) = 0;
% % %                     % a(4) = 8;
% % %                     axis(a);
% % %                     xlabel('Input SNR (dB)');
% % %                     ylabel('\gamma_0^*');
% % %
