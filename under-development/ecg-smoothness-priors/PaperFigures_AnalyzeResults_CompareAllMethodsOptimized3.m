% Similar to previous version; just changed the wavelet level calculations

close all
clc
clear

pth = 'G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES\AllResultsRawFiles\';

results_LPF = importdata([pth 'AllDataSmootherResultsNormal05_ComparisonWithLPF_Backup.txt']);
results_WAVELET = importdata([pth 'AllDataSmootherResultsNormal06_ComparisonWithWavelet2_Backup.txt']);
results_NAIVE = importdata([pth 'AllDataSmootherResultsNormal07_NaiveSweepGammaBasedOptimized_Backup.txt']);
results_FIXEDVAR = importdata([pth 'AllDataSmootherResultsNormal08_KnownVarianceOptimized_Backup.txt']);

% results_LPF = importdata([pth 'AllDataSmootherResultsArrhythmia05_ComparisonWithLPF_Backup.txt']);
% results_WAVELET = importdata([pth 'AllDataSmootherResultsArrhythmia06_ComparisonWithWavelet2_Backup.txt']);
% results_NAIVE = importdata([pth 'AllDataSmootherResultsArrhythmia07_NaiveSweepGammaBasedOptimized_Backup.txt']);
% results_FIXEDVAR = importdata([pth 'AllDataSmootherResultsArrhythmia08_KnownVarianceOptimized_Backup.txt']);

snr_list = 30 : -3: -15;

% Analyze LPF results
order_LPF = results_LPF(:, 5);
snr_LPF = results_LPF(:, 8);
snrimp_LPF = results_LPF(:, 10);
optparam_LPF = results_LPF(:, 11);

snrimp2_mean_LPF = zeros(1, length(snr_list));
snrimp2_std_LPF = zeros(1, length(snr_list));
optparam2_mean_LPF = zeros(1, length(snr_list));
optparam2_std_LPF = zeros(1, length(snr_list));

snrimp4_mean_LPF = zeros(1, length(snr_list));
snrimp4_std_LPF = zeros(1, length(snr_list));
optparam4_mean_LPF = zeros(1, length(snr_list));
optparam4_std_LPF = zeros(1, length(snr_list));

snrimp6_mean_LPF = zeros(1, length(snr_list));
snrimp6_std_LPF = zeros(1, length(snr_list));
optparam6_mean_LPF = zeros(1, length(snr_list));
optparam6_std_LPF = zeros(1, length(snr_list));
for i = 1 : length(snr_list),
    I2 = find(snr_LPF == snr_list(i) & order_LPF == 2);
    snrimp2_mean_LPF(i) = mean(snrimp_LPF(I2));
    snrimp2_std_LPF(i) = std(snrimp_LPF(I2));
    optparam2_mean_LPF(i) = mean(optparam_LPF(I2));
    optparam2_std_LPF(i) = std(optparam_LPF(I2));
    
    I4 = find(snr_LPF == snr_list(i) & order_LPF == 4);
    snrimp4_mean_LPF(i) = mean(snrimp_LPF(I4));
    snrimp4_std_LPF(i) = std(snrimp_LPF(I4));
    optparam4_mean_LPF(i) = mean(optparam_LPF(I4));
    optparam4_std_LPF(i) = std(optparam_LPF(I4));
    
    I6 = find(snr_LPF == snr_list(i) & order_LPF == 6);
    snrimp6_mean_LPF(i) = mean(snrimp_LPF(I6));
    snrimp6_std_LPF(i) = std(snrimp_LPF(I6));
    optparam6_mean_LPF(i) = mean(optparam_LPF(I6));
    optparam6_std_LPF(i) = std(optparam_LPF(I6));
end


% Analyze Naive results
snr_NAIVE = results_NAIVE(:, 8);
order_NAIVE = results_NAIVE(:, 5);
snrimp_NAIVE = results_NAIVE(:, 11);
optparam_NAIVE = results_NAIVE(:, 12);

snrimp1_mean_NAIVE = zeros(1, length(snr_list));
snrimp1_std_NAIVE = zeros(1, length(snr_list));
optparam1_mean_NAIVE = zeros(1, length(snr_list));
optparam1_std_NAIVE = zeros(1, length(snr_list));
snrimp2_mean_NAIVE = zeros(1, length(snr_list));
snrimp2_std_NAIVE = zeros(1, length(snr_list));
optparam2_mean_NAIVE = zeros(1, length(snr_list));
optparam2_std_NAIVE = zeros(1, length(snr_list));
snrimp3_mean_NAIVE = zeros(1, length(snr_list));
snrimp3_std_NAIVE = zeros(1, length(snr_list));
optparam3_mean_NAIVE = zeros(1, length(snr_list));
optparam3_std_NAIVE = zeros(1, length(snr_list));
for i = 1 : length(snr_list),
    I1 = find(snr_NAIVE == snr_list(i) & order_NAIVE == 2);
    snrimp1_mean_NAIVE(i) = mean(snrimp_NAIVE(I1));
    snrimp1_std_NAIVE(i) = std(snrimp_NAIVE(I1));
    optparam1_mean_NAIVE(i) = median(optparam_NAIVE(I1));
    optparam1_std_NAIVE(i) = std(optparam_NAIVE(I1));
    
    I2 = find(snr_NAIVE == snr_list(i) & order_NAIVE == 4);
    snrimp2_mean_NAIVE(i) = mean(snrimp_NAIVE(I2));
    snrimp2_std_NAIVE(i) = std(snrimp_NAIVE(I2));
    optparam2_mean_NAIVE(i) = median(optparam_NAIVE(I2));
    optparam2_std_NAIVE(i) = std(optparam_NAIVE(I2));
    
    I3 = find(snr_NAIVE == snr_list(i) & order_NAIVE == 6);
    snrimp3_mean_NAIVE(i) = mean(snrimp_NAIVE(I3));
    snrimp3_std_NAIVE(i) = std(snrimp_NAIVE(I3));
    optparam3_mean_NAIVE(i) = median(optparam_NAIVE(I3));
    optparam3_std_NAIVE(i) = std(optparam_NAIVE(I3));
end

% Analyze wavelet results
numeric = results_WAVELET.data;
nominal = results_WAVELET.textdata;

% TPTR = {'rigrsure'};%%%{'rigrsure', 'heursure', 'sqtwolog', 'minimaxi'};
% SORH = {'s'};%%%{'s', 'h'};
SCAL = {'sln'};%%%{'one', 'sln', 'mln'};
WLEVELS = 1 : 10;
WNAME = {'haar', 'db2', 'db3' ,'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db12', 'db16', 'coif1', 'coif2', 'coif3', 'coif4', 'coif5', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior1.5', 'bior2.6', 'bior2.8', 'bior5.5', 'bior6.8'};

snrimp_count_WAVELET = zeros(length(SCAL), length(WNAME), length(WLEVELS), length(snr_list));
snrimp_sqsum_WAVELET = zeros(length(SCAL), length(WNAME), length(WLEVELS), length(snr_list));
snrimp_sum_WAVELET = zeros(length(SCAL), length(WNAME), length(WLEVELS), length(snr_list));
wlevel_sum_WAVELET = zeros(length(SCAL), length(WNAME), length(WLEVELS), length(snr_list));
for k = 1 : length(numeric),
    scal = nominal(k, 7);
    Iscal = [];
    for j = 1 : length(SCAL),
        if(isequal(scal, SCAL(j)))
            Iscal = j;
            break;
        end
    end
    
    wname = nominal(k, 8);
    Iwname = [];
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
    
    if(~isempty(Iscal) && ~isempty(Iwname) && ~isempty(Ilevel) && ~isempty(Isnr))
        snrimp_sum_WAVELET(Iscal, Iwname, Ilevel, Isnr) = snrimp_sum_WAVELET(Iscal, Iwname, Ilevel, Isnr) + (snr1 - snr0);
        snrimp_sqsum_WAVELET(Iscal, Iwname, Ilevel, Isnr) = snrimp_sqsum_WAVELET(Iscal, Iwname, Ilevel, Isnr) + (snr1 - snr0)^2;
        wlevel_sum_WAVELET(Iscal, Iwname, Ilevel, Isnr) = wlevel_sum_WAVELET(Iscal, Iwname, Ilevel, Isnr) + level;
        snrimp_count_WAVELET(Iscal, Iwname, Ilevel, Isnr) = snrimp_count_WAVELET(Iscal, Iwname, Ilevel, Isnr) + 1;
    end
end

wlevel_mean_WAVELET = wlevel_sum_WAVELET ./ snrimp_count_WAVELET;
snrimp_mean_WAVELET = snrimp_sum_WAVELET ./ snrimp_count_WAVELET;
snrimp_std_WAVELET = sqrt(snrimp_sqsum_WAVELET./snrimp_count_WAVELET - snrimp_mean_WAVELET.^2);

snrimp_mean_averaged_WAVELET = squeeze(mean(snrimp_mean_WAVELET, 4));
snrimp_std_averaged_WAVELET = squeeze(mean(snrimp_std_WAVELET, 4));

[~, II] = sort(snrimp_mean_averaged_WAVELET(:), 'descend');
[scal_opt_WAVELET, wname_opt_WAVELET, wname_opt_WLEVELS] = ind2sub(size(snrimp_mean_averaged_WAVELET), II(1));

% Analyze Fixed-Variance Results
diff_h_FIXEDVAR = 2;%[2 4 6];
wlen_FIXEDVAR = 1; % equivalent to: [100e-3 200e-3]; % window length (s)
NoisePowerOverEstimationFactor_FIXEDVAR = [1 2 3] ; % equivalent to: [0.8 1 1.2];%%%sqrt(fs*wlen)];
ForgettingFactor_FIXEDVAR = 1;%%%0.9 : 0.1 : 1; % forgetting rate of the smoothness factor calculated in each block

wlens_FIXEDVAR = results_FIXEDVAR(:, 7);
orders_FIXEDVAR = results_FIXEDVAR(:, 8);
forgettingfactors_FIXEDVAR = results_FIXEDVAR(:, 9);
overestim_FIXEDVAR = results_FIXEDVAR(:, 10);
snr_FIXEDVAR = results_FIXEDVAR(:, 11);
snrimp1_FIXEDVAR = results_FIXEDVAR(:, 13);
snrimp2_FIXEDVAR = results_FIXEDVAR(:, 14);
optparam1_FIXEDVAR = results_FIXEDVAR(:, 15);
optparam2_FIXEDVAR = results_FIXEDVAR(:, 16);

snrimp1_mean_FIXEDVAR = zeros(length(snr_list), length(wlen_FIXEDVAR), length(diff_h_FIXEDVAR), length(ForgettingFactor_FIXEDVAR), length(NoisePowerOverEstimationFactor_FIXEDVAR));
snrimp1_std_FIXEDVAR = zeros(length(snr_list), length(wlen_FIXEDVAR), length(diff_h_FIXEDVAR), length(ForgettingFactor_FIXEDVAR), length(NoisePowerOverEstimationFactor_FIXEDVAR));
snrimp2_mean_FIXEDVAR = zeros(length(snr_list), length(wlen_FIXEDVAR), length(diff_h_FIXEDVAR), length(ForgettingFactor_FIXEDVAR), length(NoisePowerOverEstimationFactor_FIXEDVAR));
snrimp2_std_FIXEDVAR = zeros(length(snr_list), length(wlen_FIXEDVAR), length(diff_h_FIXEDVAR), length(ForgettingFactor_FIXEDVAR), length(NoisePowerOverEstimationFactor_FIXEDVAR));

% % % figure
% % % hold on;
% for i = 1 : length(snrimp1)
for ss = 1 : length(snr_list),
    for ww = 1 : length(wlen_FIXEDVAR),
        for oo = 1 : length(diff_h_FIXEDVAR),
            for ff = 1 : length(ForgettingFactor_FIXEDVAR),
                for ovr = 1 : length(NoisePowerOverEstimationFactor_FIXEDVAR),
                    I = find(snr_FIXEDVAR == snr_list(ss) & orders_FIXEDVAR == diff_h_FIXEDVAR(oo) & wlens_FIXEDVAR == wlen_FIXEDVAR(ww) & overestim_FIXEDVAR == NoisePowerOverEstimationFactor_FIXEDVAR(ovr) & forgettingfactors_FIXEDVAR == ff);
                    snrimp1_mean_FIXEDVAR(ss, ww, oo, ff, ovr) = mean(snrimp1_FIXEDVAR(I));
                    snrimp2_mean_FIXEDVAR(ss, ww, oo, ff, ovr) = mean(snrimp2_FIXEDVAR(I));
                    
                    snrimp1_std_FIXEDVAR(ss, ww, oo, ff, ovr) = std(snrimp1_FIXEDVAR(I));
                    snrimp2_std_FIXEDVAR(ss, ww, oo, ff, ovr) = std(snrimp2_FIXEDVAR(I));
                end
            end
        end
    end
end

% Plot Results
figure(1);
hold on;
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
% a = axis;
% a(1) = -9;%snr_list(end);
% a(2) = 30;%snr_list(1);
% a(3) = -3;
% a(4) = 10;
% axis(a);

% for scal = 1 : length(SCAL),
%     for wname = 1 : length(WNAME),
%         %         switch scal
%         %             case 1, col = 'r';
%         %             case 2, col = 'g';
%         %             case 3, col = 'b';
%         %         end
%         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
%         col = .6*ones(1, 3);
%         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
%         %         if(~isequal(scal, 1))
%         plot(snr_list, squeeze(snrimp_mean_WAVELET(scal, wname, :)), 'color', col);
%         %         end
%     end
% end
figure(1);
plot(snr_list, squeeze(snrimp_mean_WAVELET(scal_opt_WAVELET, wname_opt_WAVELET, wname_opt_WLEVELS, :)), 'k', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
grid

figure(2);
hold on;
set(gca, 'fontsize', 16);
set(gca, 'box', 'on');
% a = axis;
% a(1) = -9;%snr_list(end);
% a(2) = 30;%snr_list(1);
% a(3) = 0;
% a(4) = 3;
% axis(a);

% for scal = 1 : length(SCAL),
%     for wname = 1 : length(WNAME),
%         %         switch scal
%         %             case 1, col = 'r';
%         %             case 2, col = 'g';
%         %             case 3, col = 'b';
%         %         end
%         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
%         col = .6*ones(1, 3);
%         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
%         %         if(~isequal(scal, 1))
%         plot(snr_list, squeeze(snrimp_std_WAVELET(scal, wname, :)), 'color', col);
%         %         end
%     end
% end
figure(2);
plot(snr_list, squeeze(snrimp_std_WAVELET(scal_opt_WAVELET, wname_opt_WAVELET, wname_opt_WLEVELS, :)), 'k', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
grid

figure(1)
plot(snr_list, snrimp2_mean_LPF, 'c', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','c', 'MarkerFaceColor', 'c', 'MarkerSize', 6);
% plot(snr_list, snrimp4_mean_LPF, 'y', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% plot(snr_list, snrimp6_mean_LPF, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR Mean (dB)');

figure(2)
plot(snr_list, snrimp2_std_LPF, 'c', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','c', 'MarkerFaceColor', 'c', 'MarkerSize', 6);
% plot(snr_list, snrimp4_std_LPF, 'y', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% plot(snr_list, snrimp6_std_LPF, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
xlabel('Input SNR (dB)');
ylabel('\Delta SNR STD (dB)');

figure(1);
plot(snr_list, snrimp1_mean_NAIVE, 'm--', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
% plot(snr_list, snrimp2_mean_NAIVE, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
% plot(snr_list, snrimp3_mean_NAIVE, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);

figure(2);
plot(snr_list, snrimp1_std_NAIVE, 'm--', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
% plot(snr_list, snrimp2_std_NAIVE, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
% plot(snr_list, snrimp3_std_NAIVE, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);

figure(1)
for ww = 1 : length(wlen_FIXEDVAR),
    for oo = 1 : length(diff_h_FIXEDVAR),
        for ff = 1 : length(ForgettingFactor_FIXEDVAR),
            for ovr = 1 : length(NoisePowerOverEstimationFactor_FIXEDVAR),
                switch ovr
                    case 1, col = 'r';
                    case 2, col = 'g';
                    case 3, col = 'b';
                end
                %                 plot(snr_list, squeeze(snrimp1_mean_FIXEDVAR(:, ww, oo, ff, ovr)), 'm', 'marker', 'square', 'LineWidth', 1, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
                plot(snr_list, squeeze(snrimp2_mean_FIXEDVAR(:, ww, oo, ff, ovr)), 'color', col, 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor',col, 'MarkerFaceColor', col, 'MarkerSize', 4);
            end
        end
    end
end

figure(2)
for ww = 1 : length(wlen_FIXEDVAR),
    for oo = 1 : length(diff_h_FIXEDVAR),
        for ff = 1 : length(ForgettingFactor_FIXEDVAR),
            for ovr = 1 : length(NoisePowerOverEstimationFactor_FIXEDVAR),
                switch ovr
                    case 1, col = 'r';
                    case 2, col = 'g';
                    case 3, col = 'b';
                end
                %                 plot(snr_list, squeeze(snrimp1_std_FIXEDVAR(:, ww, oo, ff, ovr)), 'm', 'marker', 'square', 'LineWidth', 1, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
                plot(snr_list, squeeze(snrimp2_std_FIXEDVAR(:, ww, oo, ff, ovr)), 'color', col, 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor',col, 'MarkerFaceColor', col, 'MarkerSize', 4);
            end
        end
    end
end

snrimp1_mean_averaged_FIXEDVAR = mean(snrimp1_mean_FIXEDVAR, 1);
[~, II] = sort(snrimp1_mean_averaged_FIXEDVAR(:), 'descend');
[ww, oo, ff, ovr] = ind2sub(size(snrimp1_mean_averaged_FIXEDVAR), II(1));
% figure(1);
% plot(snr_list, squeeze(snrimp1_mean_FIXEDVAR(:, ww, oo, ff, ovr)), 'm', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% figure(2);
% plot(snr_list, squeeze(snrimp1_std_FIXEDVAR(:, ww, oo, ff, ovr)), 'm', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);

snrimp2_mean_averaged_FIXEDVAR = mean(snrimp2_mean_FIXEDVAR, 1);
[~, II] = sort(snrimp2_mean_averaged_FIXEDVAR(:), 'descend');
[ww, oo, ff, ovr] = ind2sub(size(snrimp2_mean_averaged_FIXEDVAR), II(1));
% figure(1);
% plot(snr_list, squeeze(snrimp2_mean_FIXEDVAR(:, ww, oo, ff, ovr)), 'y', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% figure(2);
% plot(snr_list, squeeze(snrimp2_std_FIXEDVAR(:, ww, oo, ff, ovr)), 'y', 'marker', 'square', 'LineWidth', 3, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);


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



% % % % % results2 = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB05_ComparisonWithLPF_Backup.txt');
% % % % % results2 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal05_ComparisonWithLPF_Backup.txt');
% % % % % results2 = [zeros(size(results.data(:,1),1), 2) results.data(:,1:end)];
% % % % 
% % % % 
% % % % % results0 = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB06_ComparisonWithWavelet2_Backup.txt'); fs = 1000;
% % % % results0 = importdata('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia06_ComparisonWithWavelet2_Backup.txt'); fs = 360;
% % % % % results0 = importdata('J:\ECGData\Normal\AllDataSmootherResultsNormal06_ComparisonWithWavelet2_Backup.txt'); fs = 128;
% % % % 
% % % % % results1 = load('J:\ECGData\PTB\AllDataSmootherResultsPTB07_NaiveSweepGammaBasedOptimized_Backup.txt');
% % % % results1 = load('J:\ECGData\Arrhythmia\AllDataSmootherResultsArrhythmia04_NaiveSweepGammaBased_Backup.txt');
% % % % % results1 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal04_NaiveSweepGammaBased_Backup.txt');
% % % % % results1 = load('J:\ECGData\Normal\AllDataSmootherResultsNormal07_NaiveSweepGammaBasedOptimized.txt');
% % % % 
% % % % 
% % % % 
% % % % numeric = results0.data(:,1:end);
% % % % nominal = results0.textdata(:,1:end);
% % % % 
% % % % snr_list = 30 : -3: -15;
% % % % % snr_list = 30;
% % % % TPTR = {'rigrsure'};%%%{'rigrsure', 'heursure', 'sqtwolog', 'minimaxi'};
% % % % SORH = {'s'};%%%{'s', 'h'};
% % % % SCAL = {'sln'};%{'one', 'sln', 'mln'};
% % % % WLEVELS = 1 : 10;
% % % % WNAME = {'haar', 'db2', 'db3' ,'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db12', 'db16', 'coif1', 'coif2', 'coif3', 'coif4', 'coif5', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior1.5', 'bior2.6', 'bior2.8', 'bior5.5', 'bior6.8'};
% % % % 
% % % % 
% % % % % % % snrmask = false(1, length(numeric));
% % % % % % % for k = 1 : length(snr_list),
% % % % % % %     snrmask(numeric(:, 3) == snr_list(k)) = true(1);
% % % % % % % end
% % % % % % %
% % % % % % % scalmask = false(1, length(numeric));
% % % % % % % for k = 1 : length(SCAL),
% % % % % % %     for j = 1 : length(nominal(:, 7)),
% % % % % % %         if(isequal(nominal(j, 7), SCAL(k)))
% % % % % % %             scalmask(j) = true(1);
% % % % % % %         end
% % % % % % %     end
% % % % % % % end
% % % % % % %
% % % % % % % numeric = numeric(snrmask & scalmask, :);
% % % % % % % nominal = nominal(snrmask & scalmask, :);
% % % % 
% % % % snrimp_count = zeros(length(SCAL), length(WNAME), length(snr_list));
% % % % snrimp_sqsum = zeros(length(SCAL), length(WNAME), length(snr_list));
% % % % snrimp_sum = zeros(length(SCAL), length(WNAME), length(snr_list));
% % % % wlevel_sum = zeros(length(SCAL), length(WNAME), length(snr_list));
% % % % for k = 1 : length(numeric),
% % % %     scal = nominal(k, 7);
% % % %     for j = 1 : length(SCAL),
% % % %         if(isequal(scal, SCAL(j)))
% % % %             Iscal = j;
% % % %             break;
% % % %         end
% % % %     end
% % % %     
% % % %     wname = nominal(k, 8);
% % % %     for j = 1 : length(WNAME),
% % % %         if(isequal(wname, WNAME(j)))
% % % %             Iwname = j;
% % % %             break;
% % % %         end
% % % %     end
% % % %     
% % % %     level = numeric(k, 1);
% % % %     Ilevel = find(WLEVELS == level, 1);
% % % %     
% % % %     snr = numeric(k, 3);
% % % %     Isnr = find(snr_list == snr, 1);
% % % %     
% % % %     snr0 = numeric(k, 4);
% % % %     snr1 = numeric(k, 5);
% % % %     
% % % %     snrimp_sum(Iscal, Iwname, Isnr) = snrimp_sum(Iscal, Iwname, Isnr) + (snr1 - snr0);
% % % %     snrimp_sqsum(Iscal, Iwname, Isnr) = snrimp_sqsum(Iscal, Iwname, Isnr) + (snr1 - snr0)^2;
% % % %     wlevel_sum(Iscal, Iwname, Isnr) = wlevel_sum(Iscal, Iwname, Isnr) + level;
% % % %     snrimp_count(Iscal, Iwname, Isnr) = snrimp_count(Iscal, Iwname, Isnr) + 1;
% % % % end
% % % % 
% % % % wlevel_mean = wlevel_sum ./ snrimp_count;
% % % % snrimp_mean = snrimp_sum ./ snrimp_count;
% % % % snrimp_std = sqrt(snrimp_sqsum./snrimp_count - snrimp_mean.^2);
% % % % 
% % % % snrimp_mean_averaged = squeeze(mean(snrimp_mean, 3));
% % % % snrimp_std_averaged = squeeze(mean(snrimp_std, 3));
% % % % 
% % % % [YY II] = sort(snrimp_mean_averaged(:), 'descend');
% % % % 
% % % % [scal_opt, wname_opt] = ind2sub(size(snrimp_mean_averaged), II(1));
% % % % 
% % % % % [YY II] = sort(snrimp_std_averaged(:), 'ascend');
% % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % % % curve_name = cell(length(SCAL), length(WNAME));
% % % % % % % for scal = 1 : length(SCAL),
% % % % % % %     for wname = 1 : length(WNAME),
% % % % % % %         curve_name(scal, wname) = {[SCAL{scal} '_' WNAME{wname} ' ' num2str(snrimp_mean_averaged(scal,wname)) '(' num2str(snrimp_std_averaged(scal,wname)) ')']};
% % % % % % %         switch scal
% % % % % % %             case 1, col = 'r';
% % % % % % %             case 2, col = 'g';
% % % % % % %             case 3, col = 'b';
% % % % % % %         end
% % % % % % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % % % % % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % % % % % %         plot(snr_list, squeeze(wlevel_mean(scal, wname, :)), 'color', col);
% % % % % % %     end
% % % % % % % end
% % % % % % % grid
% % % % % % % curve_name(II)
% % % % 
% % % % figure
% % % % hold on
% % % % for scal = 1 : length(SCAL),
% % % %     for wname = 1 : length(WNAME),
% % % %         %         switch scal
% % % %         %             case 1, col = 'r';
% % % %         %             case 2, col = 'g';
% % % %         %             case 3, col = 'b';
% % % %         %         end
% % % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % % %         col = .6*ones(1, 3);
% % % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % % %         %         if(~isequal(scal, 1))
% % % %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), 'color', col);
% % % %         %         end
% % % %     end
% % % % end
% % % % plot(snr_list, squeeze(snrimp_mean(scal_opt, wname_opt, :)), 'k', 'linewidth', 3);
% % % % % grid
% % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % % % for scal = 1 : length(SCAL),
% % % % % % %     for wname = 1 : length(WNAME),
% % % % % % %         switch scal
% % % % % % %             case 1, col = 'r';
% % % % % % %             case 2, col = 'g';
% % % % % % %             case 3, col = 'b';
% % % % % % %         end
% % % % % % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % % % % % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % % % % % %         plot(snr_list, squeeze(snrimp_std(scal, wname, :)), 'color', col);
% % % % % % %     end
% % % % % % % end
% % % % % % % grid
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % for scal = 1 : length(SCAL),
% % % % % % %     for wname = 1 : length(WNAME),
% % % % % % %         switch scal
% % % % % % %             case 1, col = 'r';
% % % % % % %             case 2, col = 'g';
% % % % % % %             case 3, col = 'b';
% % % % % % %         end
% % % % % % %         %         col = 0.6*scal*ones(1, 3)/length(SCAL);
% % % % % % %         %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
% % % % % % %         errorbar(snr_list, squeeze(snrimp_mean(scal, wname, :)), squeeze(snrimp_std(scal, wname, :)), 'color', col);
% % % % % % %     end
% % % % % % % end
% % % % % % % grid
% % % % 
% % % % 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 
% % % % snr_list = 30 : -3: -15;
% % % % 
% % % % snr = results1(:, 8);
% % % % order = results1(:, 5);
% % % % snr0 = results1(:, 9);
% % % % % snrimp1 = results1(:, 10);
% % % % % optparam = results1(:, 11);
% % % % snrimp1 = results1(:, 11);
% % % % optparam = results1(:, 12);
% % % % 
% % % % snrimp1_mean = zeros(1, length(snr_list));
% % % % snrimp1_std = zeros(1, length(snr_list));
% % % % optparam1_mean = zeros(1, length(snr_list));
% % % % optparam1_std = zeros(1, length(snr_list));
% % % % snrimp2_mean = zeros(1, length(snr_list));
% % % % snrimp2_std = zeros(1, length(snr_list));
% % % % optparam2_mean = zeros(1, length(snr_list));
% % % % optparam2_std = zeros(1, length(snr_list));
% % % % snrimp3_mean = zeros(1, length(snr_list));
% % % % snrimp3_std = zeros(1, length(snr_list));
% % % % optparam3_mean = zeros(1, length(snr_list));
% % % % optparam3_std = zeros(1, length(snr_list));
% % % % for i = 1 : length(snr_list),
% % % %     I1 = find(snr == snr_list(i) & order == 2);
% % % %     snrimp1_mean(i) = mean(snrimp1(I1));
% % % %     snrimp1_std(i) = std(snrimp1(I1));
% % % %     optparam1_mean(i) = median(optparam(I1));
% % % %     optparam1_std(i) = std(optparam(I1));
% % % %     
% % % %     I2 = find(snr == snr_list(i) & order == 4);
% % % %     snrimp2_mean(i) = mean(snrimp1(I2));
% % % %     snrimp2_std(i) = std(snrimp1(I2));
% % % %     optparam2_mean(i) = median(optparam(I2));
% % % %     optparam2_std(i) = std(optparam(I2));
% % % %     
% % % %     I3 = find(snr == snr_list(i) & order == 6);
% % % %     snrimp3_mean(i) = mean(snrimp1(I3));
% % % %     snrimp3_std(i) = std(snrimp1(I3));
% % % %     optparam3_mean(i) = median(optparam(I3));
% % % %     optparam3_std(i) = std(optparam(I3));
% % % % end
% % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% % % % % errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% % % % % plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% % % % % plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
% % % % plot(snr_list, snrimp1_mean, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','b', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
% % % % plot(snr_list, snrimp2_mean, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
% % % % plot(snr_list, snrimp3_mean, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
% % % % % plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % grid
% % % % % % % set(gca, 'fontsize', 16);
% % % % % % % set(gca, 'box', 'on');
% % % % % % % axis tight
% % % % % % % a = axis;
% % % % % % % a(1) = -16;%snr_list(end);
% % % % % % % a(2) = 31;%snr_list(1);
% % % % % % % % a(3) = 0;
% % % % % % % % a(4) = 8;
% % % % % % % axis(a);
% % % % % % % xlabel('Input SNR (dB)');
% % % % % % % ylabel('\Delta SNR Mean (dB)');
% % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % % % plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % % % % plot(snr_list, snrimp3_std, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % % % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % grid
% % % % % % % set(gca, 'fontsize', 16);
% % % % % % % set(gca, 'box', 'on');
% % % % % % % set(gcf, 'Position', [520 378 560 219])
% % % % % % % axis tight
% % % % % % % a = axis;
% % % % % % % a(1) = -16;%snr_list(end);
% % % % % % % a(2) = 31;%snr_list(1);
% % % % % % % % a(3) = 0;
% % % % % % % % a(4) = 8;
% % % % % % % axis(a);
% % % % % % % xlabel('Input SNR (dB)');
% % % % % % % ylabel('\Delta SNR SD (dB)');
% % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % % % % errorbar(snr_list, optparam1_mean, optparam1_std, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % % errorbar(snr_list, optparam2_mean, optparam2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % % errorbar(snr_list, optparam3_mean, optparam3_std, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % plot(snr_list, optparam1_mean, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % % % % % plot(snr_list, optparam2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % % % % % plot(snr_list, optparam3_mean, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % % % % % set(gca, 'yscale', 'log');
% % % % % % % grid
% % % % % % % set(gca, 'fontsize', 16);
% % % % % % % set(gca, 'box', 'on');
% % % % % % % set(gcf, 'Position', [520 81 560 219])
% % % % % % % axis tight
% % % % % % % a = axis;
% % % % % % % a(1) = -16;%snr_list(end);
% % % % % % % a(2) = 31;%snr_list(1);
% % % % % % % % a(3) = 0;
% % % % % % % % a(4) = 8;
% % % % % % % axis(a);
% % % % % % % xlabel('Input SNR (dB)');
% % % % % % % ylabel('\gamma_0^*');
% % % % 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 
% % % % % snr_list = 30 : -3: -12;
% % % % snr_list = 30 : -3: -15;
% % % % 
% % % % snr = results2(:, 8);
% % % % snr0 = results2(:, 9);
% % % % snrimp1 = results2(:, 10);
% % % % % snrimp2 = results(:, 10);
% % % % % snrimp3 = results(:, 11);
% % % % optparam = results2(:, 11);
% % % % 
% % % % snrimp1_mean = zeros(1, length(snr_list));
% % % % % snrimp2_mean = zeros(1, length(snr_list));
% % % % % snrimp3_mean = zeros(1, length(snr_list));
% % % % snrimp1_std = zeros(1, length(snr_list));
% % % % % snrimp2_std = zeros(1, length(snr_list));
% % % % % snrimp3_std = zeros(1, length(snr_list));
% % % % optparam_mean = zeros(1, length(snr_list));
% % % % optparam_std = zeros(1, length(snr_list));
% % % % for i = 1 : length(snr_list),
% % % %     I = find(snr == snr_list(i));
% % % %     snrimp1_mean(i) = mean(snrimp1(I));
% % % %     %     snrimp2_mean(i) = mean(snrimp2(I));
% % % %     %     snrimp3_mean(i) = mean(snrimp3(I));
% % % %     snrimp1_std(i) = std(snrimp1(I));
% % % %     %     snrimp2_std(i) = std(snrimp2(I));
% % % %     %     snrimp3_std(i) = std(snrimp3(I));
% % % %     optparam_mean(i) = mean(optparam(I));
% % % %     optparam_std(i) = std(optparam(I));
% % % % end
% % % % 
% % % % % figure
% % % % % hold on
% % % % % errorbar(snr_list, snrimp1_mean, snrimp1_std, 'b', 'linewidth', 2.5);
% % % % % errorbar(snr_list, snrimp2_mean, snrimp2_std, 'r', 'linewidth', 2);
% % % % % plot(snr_list, snrimp1_mean, 'b', 'linewidth', 2.5);
% % % % % plot(snr_list, snrimp2_mean, 'r', 'linewidth', 2);
% % % % plot(snr_list, snrimp1_mean, 'c', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','m', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% % % % % plot(snr_list, snrimp2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % grid
% % % % % % % set(gca, 'fontsize', 16);
% % % % % % % set(gca, 'box', 'on');
% % % % % % % a = axis;
% % % % % % % a(1) = -15;%snr_list(end);
% % % % % % % a(2) = 33;%snr_list(1);
% % % % % % % a(3) = 0;
% % % % % % % a(4) = 8;
% % % % % axis(a);
% % % % xlabel('Input SNR (dB)');
% % % % ylabel('\Delta SNR Mean (dB)');
% % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % % % plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 6);
% % % % % % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % grid
% % % % % % % set(gca, 'fontsize', 16);
% % % % % % % set(gca, 'box', 'on');
% % % % % % % a = axis;
% % % % % % % a(1) = -15;%snr_list(end);
% % % % % % % a(2) = 33;%snr_list(1);
% % % % % % % a(3) = 0;
% % % % % % % a(4) = 2;
% % % % % % % % axis(a);
% % % % % % % xlabel('Input SNR (dB)');
% % % % % % % ylabel('\Delta SNR SD (dB)');
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % errorbar(snr_list, optparam_mean, optparam_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % grid
% % % % % % % set(gca, 'fontsize', 16);
% % % % % % % set(gca, 'box', 'on');
% % % % % % % % a = axis;
% % % % % % % % a(1) = -15;%snr_list(end);
% % % % % % % % a(2) = 33;%snr_list(1);
% % % % % % % % a(3) = 0;
% % % % % % % % a(4) = 2;
% % % % % % % % axis(a);
% % % % % % % xlabel('Input SNR (dB)');
% % % % % % % ylabel('Optim Param');
% % % % % % %

% % % % % % % figure
% % % % % % % hold on
% % % % % % % plot(snr_list, snrimp1_std, 'b', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % % % % plot(snr_list, snrimp3_std, 'g', 'marker', 'square', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
% % % % % % % % plot(snr_list, snrimp2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % grid
% % % % % % % 
% % % % % % % figure
% % % % % % % hold on
% % % % % % % % errorbar(snr_list, optparam1_mean, optparam1_std, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % % errorbar(snr_list, optparam2_mean, optparam2_std, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % % errorbar(snr_list, optparam3_mean, optparam3_std, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 6);
% % % % % % % plot(snr_list, optparam1_mean, 'b', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % % % % % plot(snr_list, optparam2_mean, 'r', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % % % % % plot(snr_list, optparam3_mean, 'g', 'marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize', 5);
% % % % % % % set(gca, 'yscale', 'log');
% % % % % % % grid

