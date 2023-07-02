close all
clc
clear
% results = importdata('J:\ECGData\PTB\AllDataSmootherResultsPTB06_ComparisonWithWavelet2_Backup.txt'); fs = 1000;
% results = importdata('G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES\AllResultsRawFiles\AllDataSmootherResultsArrhythmia06_ComparisonWithWavelet2_Backup.txt'); fs = 360;
results = importdata('G:\Sameni\Documents\Papers\Journal 40 (Spline ECG Filtering)\CODES\AllResultsRawFiles\AllDataSmootherResultsNormal06_ComparisonWithWavelet2_Backup.txt'); fs = 128;
numeric = results.data(:,1:end);
nominal = results.textdata(:,1:end);

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

snrimp_mean_averaged = mean(snrimp_mean, 3);
snrimp_std_averaged = mean(snrimp_std, 3);

[YY II] = sort(snrimp_mean_averaged(:), 'descend');
% [YY II] = sort(snrimp_std_averaged(:), 'ascend');

figure
hold on
curve_name = cell(length(SCAL), length(WNAME));
for scal = 1 : length(SCAL),
    for wname = 1 : length(WNAME),
        curve_name(scal, wname) = {[SCAL{scal} '_' WNAME{wname} ' ' num2str(snrimp_mean_averaged(scal,wname)) '(' num2str(snrimp_std_averaged(scal,wname)) ')']};
        switch scal
            case 1, col = 'r';
            case 2, col = 'g';
            case 3, col = 'b';
        end
        %         col = 0.6*scal*ones(1, 3)/length(SCAL);
        %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
        plot(snr_list, squeeze(wlevel_mean(scal, wname, :)), 'color', col);
    end
end
grid
curve_name(II)

figure
hold on
for scal = 1 : length(SCAL),
    for wname = 1 : length(WNAME),
        switch scal
            case 1, col = 'r';
            case 2, col = 'g';
            case 3, col = 'b';
        end
        %         col = 0.6*scal*ones(1, 3)/length(SCAL);
        %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
        plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), 'color', col);
    end
end
grid

figure
hold on
for scal = 1 : length(SCAL),
    for wname = 1 : length(WNAME),
        switch scal
            case 1, col = 'r';
            case 2, col = 'g';
            case 3, col = 'b';
        end
        %         col = 0.6*scal*ones(1, 3)/length(SCAL);
        %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
        plot(snr_list, squeeze(snrimp_std(scal, wname, :)), 'color', col);
    end
end
grid

figure
hold on
for scal = 1 : length(SCAL),
    for wname = 1 : length(WNAME),
        switch scal
            case 1, col = 'r';
            case 2, col = 'g';
            case 3, col = 'b';
        end
        %         col = 0.6*scal*ones(1, 3)/length(SCAL);
        %         plot(snr_list, squeeze(snrimp_mean(scal, wname, :)), '.', 'color', col);
        errorbar(snr_list, squeeze(snrimp_mean(scal, wname, :)), squeeze(snrimp_std(scal, wname, :)), 'color', col);
    end
end
grid

