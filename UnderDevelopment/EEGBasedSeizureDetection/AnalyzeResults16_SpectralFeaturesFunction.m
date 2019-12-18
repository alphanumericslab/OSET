function AnalyzeResults16_SpectralFeaturesFunction(Rdata, Rtextdata, outputfileame1, outputfileame2, per_subject, ntrees)
% close all;
% % parameters
% % % % logfileame = 'SVMOptimizedParamters03.txt';
% outputfileame1 = 'Submission55.csv';
% outputfileame2 = 'Submission55_unbiased.csv';
% % % bins = [50 75 100 120 180 250 300 500 750];%75:6:81;%120;
% % % bins = 120;%50:1:550;%:13:250;%75:7:120;%75:6:81;%120;
% % % bins = 77;
% gmorder = 5;
% dropfeatures = 5;
% sigma = 150*ones(1, 7);
% sigma = [181.5 27.0 179 66 111 26 19];
% per_subject = 0;
% cross_validation_number = 5;
% train_ratio_per_cross_validation = 0.5; % the ratio of the interictal/preictal data used for test for each cross-validation iteration
% minAUCThreshold = 0.53;
% % meanAUCThreshold = 0.6;
% % maxAUCThreshold = 0.75;
% minPreictalPercentage = 40;
% maxPreictalPercentage = 60;

% % fix the random seed for randperm
% rng(0);

% load data
% % % R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
% % % R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
% % % R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
% % % R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
% % % R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');
% % % R6 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures3.txt');

filenames_all = Rtextdata;
subjects_all = Rdata(:, 1);
% trials_all = R1.data(:, 2);
types_all = Rdata(:, 3);
interictal_label = 1;
preictal_label = 2;
test_label = 3;
% dropfeatures = 0;
% r_raw = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)]) R6.data(:, 4:end)];
r_raw = Rdata(:, 4:end);

% discard/keep some features
% discardcols = [3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results
% r_raw(:, discardcols) = [];
% goodfeatures = [4 5 8 10 11 13 14 15 16 17 19 20 21 22 23 26 31 36 37 38];
% r_raw = r_raw(:, goodfeatures);

% log file
outfid1 = fopen(outputfileame1, 'w');
fprintf(outfid1,'clip,preictal\n');
outfid2 = fopen(outputfileame2, 'w');
fprintf(outfid2,'clip,preictal\n');

% per subject feature normalization
for s = 1:max(subjects_all),
    sb = (subjects_all == s);
    r_raw(sb, :) = (r_raw(sb, :) - ones(size(r_raw(sb, :), 1), 1)*mean(r_raw(sb, :), 1))./(ones(size(r_raw(sb, :), 1), 1)*std(r_raw(sb, :), [], 1)); % remove mean and normalize
    % r_raw(sb, :) = r_raw(sb, :) - ones(size(r_raw(sb, :), 1), 1)*mean(r_raw(sb, :), 1); % remove the mean only
end

% % PCA: decorrelate and reduce the features
% Cx = cov(r_raw);
% [V, D] = eig(Cx);
% D = diag(D);
% [~, I] = sort(D, 1, 'descend');
% w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues
% r_raw = r_raw*w;

% figure
% stem(D(I)/D(I(1)));
% grid

if(per_subject == 1)
    max_subject = max(subjects_all);
else
    max_subject = 1;
end

% scores = [];
for s = 1:max_subject,
    if(per_subject == 1)
        subjects = (subjects_all == s);
    else
        subjects = true(size(subjects_all));
    end
    
    filenames = filenames_all(subjects); % file names per subject
    r = r_raw(subjects, :); % all features per subject
    %     trials = trials_all(subjects); % trials per subject
    types = types_all(subjects); % types (interictal/preictal/test) per subject
    interictal_indexes = find(types == interictal_label); % interictal
    preictal_indexes = find(types == preictal_label); % preictal
    train_indexes = [interictal_indexes ; preictal_indexes];
    test_indexes = find(types == test_label); % test
    %testfiles = filenames(test_indexes);
    
    % % % % %     [estimated_class error] = classify(r, r(train_indexes, :), types(train_indexes), 'quadratic');
    b = TreeBagger(ntrees, r(train_indexes, :), types(train_indexes),'OOBPred','on','Method','classification');%,'prior', [0.5 0.5]);%, 'Cost',[0 1.0 ; 1.0 0]);
    [~, scores, ~] = predict(b, r);
    
    % % %     figure
    % % %     subplot(121);
    % % %     stem(strcmp(estimated_class(train_indexes), num2str(preictal_label)) - 0.5);
    % % %     grid
    % % %
    % % %     subplot(122);
    % % %     stem(strcmp(estimated_class(test_indexes), num2str(preictal_label)) - 0.5);
    % % %     grid
    
    
    
    
    
    subject_bias = ( median(scores(interictal_indexes, preictal_label)) + median(scores(preictal_indexes, preictal_label)) ) / 2;
    train_bias = median(scores(test_indexes, preictal_label));
    disp(['Subject Bias = ' num2str(subject_bias)]);
    disp(['median(scores(interictal_indexes, preictal_label))(std) = ' num2str(median(scores(interictal_indexes, preictal_label))) '(' num2str(std(scores(interictal_indexes, preictal_label))) ')']);
    disp(['median(scores(preictal_indexes, preictal_label))(std) = ' num2str(median(scores(preictal_indexes, preictal_label))) '(' num2str(std(scores(preictal_indexes, preictal_label))) ')']);
    disp(['median(scores(test_indexes, preictal_label))(std) = ' num2str(median(scores(test_indexes, preictal_label))) '(' num2str(std(scores(test_indexes, preictal_label))) ')']);
    for i = 1:length(test_indexes),
        % %         if(isempty(scores))
        % %             fprintf(outfid,'%s,%10.8f\n', filenames{test_indexes(i)}, double(estimated_class(test_indexes(i)) == preictal_label));
        % %         else
        % % % % % % % % %         fprintf(outfid2,'%s,%10.8f\n', filenames{test_indexes(i)}, strcmp(estimated_class(train_indexes(i)), num2str(preictal_label)));
        fprintf(outfid1,'%s,%10.8f\n', filenames{test_indexes(i)}, scores(test_indexes(i), preictal_label));
%         fprintf(outfid2,'%s,%10.8f\n', filenames{test_indexes(i)}, max(min(scores(test_indexes(i), preictal_label) - subject_bias + 0.5, 1.0), 0.0));
        fprintf(outfid2,'%s,%10.8f\n', filenames{test_indexes(i)}, max(min(scores(test_indexes(i), preictal_label) - train_bias + subject_bias, 1.0), 0.0));
        % %         end
    end
end
fclose(outfid1);
fclose(outfid2);