close all
clc
clear
% results0 = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB06_ComparisonWithWavelet2_Backup.txt'); fs = 1000;
% results0 = importdata('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia06_ComparisonWithWavelet2_Backup.txt'); fs = 360;
results0 = importdata('J:\ECGData\Normal\AllDataSmootherResultsNormal06_ComparisonWithWavelet2_Backup.txt'); fs = 128;

% results1 = load('J:\ECGData\PTB\AllDataSmootherResultsPTB04_NaiveSweepGammaBased_Backup.txt');
% results1 = load('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia04_NaiveSweepGammaBased_Backup.txt');
results1 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal04_NaiveSweepGammaBased_Backup.txt');

% results2 = load('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia05_ComparisonWithLPF_Backup.txt');
% results2 = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB05_ComparisonWithLPF_Backup.txt');
results2 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal05_ComparisonWithLPF_Backup.txt');
% results2 = [zeros(size(results.data(:,1),1), 2) results.data(:,1:end)];



numeric = results0.data(:,1:end);
nominal = results0.textdata(:,1:end);

snr_list = 30 : -3: -15;
% snr_list = 30;
TPTR = {'rigrsure'};%%%{'rigrsure', 'heursure', 'sqtwolog', 'minimaxi'};
SORH = {'s'};%%%{'s', 'h'};
SCAL = {'one', 'sln', 'mln'};
WLEVELS = 1 : 10;
WNAME = {'haar', 'db2', 'db3' ,'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db12', 'db16', 'coif1', 'coif2', 'coif3', 'coif4', 'coif5', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior1.5', 'bior2.6', 'bior2.8', 'bior5.5', 'bior6.8'};


% % % snrmask = false(1, length(numeric));
% % % for k = 1 : length(snr_list),
% % %     snrmask(numeric(:, 3) == snr_list(k)) = true(1);
% % % end
% % %
% % % scalmask = false(1, length(numeric));
% % % for k = 1 : length(SCAL),
% % %     for j = 1 : length(nominal(:, 7)),
% % %         if(isequal(nominal(j, 7), SCAL(k)))
% % %             scalmask(j) = true(1);
% % %         end
% % %     end
% % % end
% % %
% % % numeric = numeric(snrmask & scalmask, :);
% % % nominal = nominal(snrmask & scalmask, :);

snrimp_count = zeros(length(SCAL), length(WNAME), length(snr_list));
snrimp_sqsum = zeros(length(SCAL), length(WNAME), length(snr_list));
snrimp_sum = zeros(length(SCAL), length(WNAME), length(snr_list));
wlevel_sum = zeros(length(SCAL), length(WNAME), length(snr_list));
for k = 1 : length(numeric),
    scal = nominal(k, 7);
    for j = 1 : length(SCAL),
        if(isequal(scal, SCAL(j)))
            Iscal = j;
            break;
        end
    end
    
    wname = nominal(k, 8);
    for j = 1 : length(WNAME),
        if(isequal(wname, WNAME(j)))
            Iwname = j;
            break;
        end
    end
    
    level = numeric(k, 1);
    Ilevel = find(WLEVELS == level, 1);
    
    snr = numeric(k, 3);
    Isnr = find(snr_list == snr, 1);
    
    snr0 = numeric(k, 4);
    snr1 = numeric(k, 5);
    
    snrimp_sum(Iscal, Iwname, Isnr) = snrimp_sum(Iscal, Iwname, Isnr) + (snr1 - snr0);
    snrimp_sqsum(Iscal, Iwname, Isnr) = snrimp_sqsum(Iscal, Iwname, Isnr) + (snr1 - snr0)^2;
    wlevel_sum(Iscal, Iwname, Isnr) = wlevel_sum(Iscal, Iwname, Isnr) + level;
    snrimp_count(Iscal, Iwname, Isnr) = snrimp_count(Iscal, Iwname, Isnr) + 1;
end

wlevel_mean = wlevel_sum ./ snrimp_count;
snrimp_mean = snrimp_sum ./ snrimp_count;
snrimp_std = sqrt(snrimp_sqsum./snrimp_count - snrimp_mean.^2);

snrimp_mean_averaged = squeeze(mean(snrimp_mean, 3));
snrimp_std_averaged = squeeze(mean(snrimp_std, 3));

[YY II] = sort(snrimp_mean_averaged(:), 'descend');

[scal_opt, wname_opt] = ind2sub(size(snrimp_mean_averaged), II(1));

% [YY II] = sort(snrimp_std_averaged(:), 'ascend');

% % % figure
% % % hold on
% % % curve_name = cell(length(SCAL), length(WNAME));
% % % for scal = 1 : length(SCAL),
% % %     for wname = 1 : length(WNAME),
% % %         curve_name(scal, wname) = {[SCAL{scal} '_' WNAME{wname} ' ' num2str(snrimp_mean_averaged(scal,wname)) '(' num2str(snrimp_std_averaged(scal,wname)) ')']};
% % %         switch scal
% % %             case 1, col = 'r';
% % %             case 2, col = 'g';
% % %             case 3, col = 'b';
% % %         end
% % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % %         plot(snr_list, squeeze(wlevel_mean(scal, wname, :)), 'color', col);
% % %     end
% % % end
% % % grid
% % % curve_name(II)

figure
hold on
for scal = 1 : length(SCAL),
    for wname = 1 : length(WNAME),
        %         switch scal
        %             case 1, col = 'r';
        %             case 2, col = 'g';
        %             case 3, col = 'b';
        %         end
        %         col = 0.6*scal*ones(1, 3)/length(SCAL);
        col = .6*ones(1, 3);
        %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
        %         if(~isequal(scal, 1))
        plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), 'color', col);
        %         end
    end
end
plot(snr_list, squeeze(snrimp_mean(scal_opt, wname_opt, :)), 'k', 'linewidth', 3);
% grid

% % % figure
% % % hold on
% % % for scal = 1 : length(SCAL),
% % %     for wname = 1 : length(WNAME),
% % %         switch scal
% % %             case 1, col = 'r';
% % %             case 2, col = 'g';
% % %             case 3, col = 'b';
% % %         end
% % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % %         plot(snr_list, squeeze(snrimp_std(scal, wname, :)), 'color', col);
% % %     end
% % % end
% % % grid
% % %
% % % figure
% % % hold on
% % % for scal = 1 : length(SCAL),
% % %     for wname = 1 : length(WNAME),
% % %         switch scal
% % %             case 1, col = 'r';
% % %             case 2, col = 'g';
% % %             case 3, col = 'b';
% % %         end
% % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % %         errorbar(snr_list, squeeze(snrimp_mean(scal, wname, :)), squeeze(snrimp_std(scal, wname, :)), 'color', col);
% % %     end
% % % end
% % % grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr_list = 30 : -3: -15;

snr = results1(:, 8);
order = results1(:, 5);
snr0 = results1(:, 9);
snrimp1 = results1(:, 10);
optparam = results1(:, 11);

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

% % % figure
% % % hold on
% errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
plot(snr_list, snrimp1_mean, 'b', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','b', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
plot(snr_list, snrimp2_mean, 'r', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
plot(snr_list, snrimp3_mean, 'g', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
% plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
grid
% % % set(gca, 'fontsize', 16);
% % % set(gca, 'box', 'on');
% % % axis tight
% % % a = axis;
% % % a(1) = -16;%snr_list(end);
% % % a(2) = 31;%snr_list(1);
% % % % a(3) = 0;
% % % % a(4) = 8;
% % % axis(a);
% % % xlabel('Input SNR (dB)');
% % % ylabel('\Delta SNR Mean (dB)');

% % % figure
% % % hold on
% % % plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % plot(snr_list, snrimp2_std, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % plot(snr_list, snrimp3_std, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % grid
% % % set(gca, 'fontsize', 16);
% % % set(gca, 'box', 'on');
% % % set(gcf, 'Position', [520 378 560 219])
% % % axis tight
% % % a = axis;
% % % a(1) = -16;%snr_list(end);
% % % a(2) = 31;%snr_list(1);
% % % % a(3) = 0;
% % % % a(4) = 8;
% % % axis(a);
% % % xlabel('Input SNR (dB)');
% % % ylabel('\Delta SNR SD (dB)');

% % % figure
% % % hold on
% % % % errorbar(snr_list, optparam1_mean, optparam1_std, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % errorbar(snr_list, optparam2_mean, optparam2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % errorbar(snr_list, optparam3_mean, optparam3_std, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % plot(snr_list, optparam1_mean, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % plot(snr_list, optparam2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % plot(snr_list, optparam3_mean, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % set(gca, 'yscale', 'log');
% % % grid
% % % set(gca, 'fontsize', 16);
% % % set(gca, 'box', 'on');
% % % set(gcf, 'Position', [520 81 560 219])
% % % axis tight
% % % a = axis;
% % % a(1) = -16;%snr_list(end);
% % % a(2) = 31;%snr_list(1);
% % % % a(3) = 0;
% % % % a(4) = 8;
% % % axis(a);
% % % xlabel('Input SNR (dB)');
% % % ylabel('\gamma_0^*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% snr_list = 30 : -3: -12;
snr_list = 30 : -3: -15;

snr = results2(:, 8);
snr0 = results2(:, 9);
snrimp1 = results2(:, 10);
% snrimp2 = results(:, 10);
% snrimp3 = results(:, 11);
optparam = results2(:, 11);

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

% figure
% hold on
% errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
plot(snr_list, snrimp1_mean, 'c', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % grid
% % % set(gca, 'fontsize', 16);
% % % set(gca, 'box', 'on');
% % % a = axis;
% % % a(1) = -15;%snr_list(end);
% % % a(2) = 33;%snr_list(1);
% % % a(3) = 0;
% % % a(4) = 8;
% axis(a);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR Mean (dB)');

% % % figure
% % % hold on
% % % plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % grid
% % % set(gca, 'fontsize', 16);
% % % set(gca, 'box', 'on');
% % % a = axis;
% % % a(1) = -15;%snr_list(end);
% % % a(2) = 33;%snr_list(1);
% % % a(3) = 0;
% % % a(4) = 2;
% % % % axis(a);
% % % xlabel('Input SNR (dB)');
% % % ylabel('\Delta SNR SD (dB)');
% % %
% % % figure
% % % hold on
% % % errorbar(snr_list, optparam_mean, optparam_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % grid
% % % set(gca, 'fontsize', 16);
% % % set(gca, 'box', 'on');
% % % % a = axis;
% % % % a(1) = -15;%snr_list(end);
% % % % a(2) = 33;%snr_list(1);
% % % % a(3) = 0;
% % % % a(4) = 2;
% % % % axis(a);
% % % xlabel('Input SNR (dB)');
% % % ylabel('Optim Param');
% % %
