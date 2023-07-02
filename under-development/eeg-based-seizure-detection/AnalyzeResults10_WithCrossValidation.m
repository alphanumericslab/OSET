clear
close all;
clc

% parameters
ofname = 'SVMOptimizedParamters02.txt';
% % % bins = [50 75 100 120 180 250 300 500 750];%75:6:81;%120;
% % % bins = 120;%50:1:550;%:13:250;%75:7:120;%75:6:81;%120;
% % % bins = 77;
% gmorder = 5;
% dropfeatures = 5;
sigma = 150*ones(1, 7);
% sigma = [181.5 27.0 179 66 111 26 19];
per_subject = 1;
cross_validation_number = 10;
train_ratio_per_cross_validation = 0.5; % the ratio of the interictal/preictal data used for test for each cross-validation iteration

% fix the random seed for randperm
rng(0);

% load data
R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');

filenames_all = R1.textdata;
subjects_all = R1.data(:, 1);
trials_all = R1.data(:, 2);
types_all = R1.data(:, 3);
interictal_label = 1;
preictal_label = 2;
test_label = 3;
r_raw = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)])];

% discard/keep some features
% discardcols = [3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results
% r_raw(:, discardcols) = [];
% goodfeatures = [4 5 8 10 11 13 14 15 16 17 19 20 21 22 23 26 31 36 37 38];
% r_raw = r_raw(:, goodfeatures);

% log file
fid = fopen(ofname,'w');

% per subject feature normalization
for s = 1:max(subjects_all),
    sb = (subjects_all == s);
    r_raw(sb, :) = (r_raw(sb, :) - ones(size(r_raw(sb, :), 1), 1)*mean(r_raw(sb, :), 1))./(ones(size(r_raw(sb, :), 1), 1)*std(r_raw(sb, :), [], 1)); % remove mean and normalize
    % r_raw(sb, :) = r_raw(sb, :) - ones(size(r_raw(sb, :), 1), 1)*mean(r_raw(sb, :), 1); % remove the mean only
end

% % % % PCA: decorrelate and reduce the features
% % % Cx = cov(r_raw);
% % % [V, D] = eig(Cx);
% % % D = diag(D);
% % % [~, I] = sort(D, 1, 'descend');
% % % w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues
% % % r_raw = r_raw*w;

if(per_subject == 1)
    max_subject = max(subjects_all);
else
    max_subject = 1;
end

for s = 1:max_subject,
    if(per_subject == 1)
        subjects = (subjects_all == s);
    else
        subjects = true(size(subjects_all));
    end
    
    filenames = filenames_all(subjects); % file names per subject
    r = r_raw(subjects, :); % all features per subject
    trials = trials_all(subjects); % trials per subject
    types = types_all(subjects); % types (interictal/preictal/test) per subject
    interictal_indexes = find(types == interictal_label); % interictal
    preictal_indexes = find(types == preictal_label); % preictal
    test_indexes = find(types == test_label); % test
    
    % cross validate over labeled data
    for f = [1./(10:-1:1) (2:10)],%[1/10 1/9 1/8 1/5 1/4 1/3 1/2 1 2 3 4 5 10],
        for p = -149:1:149,
            AUC = zeros(cross_validation_number, 1);
            PreictalPercentage = zeros(cross_validation_number, 1);
            for k = 1:cross_validation_number,
                prm_i = randperm(length(interictal_indexes));
                tr_i = floor(train_ratio_per_cross_validation*length(interictal_indexes));
                interictal_train_index_subset = interictal_indexes(prm_i(1:tr_i));
                
                prm_p = randperm(length(preictal_indexes));
                tr_p = floor(train_ratio_per_cross_validation*length(preictal_indexes));
                preictal_train_index_subset = preictal_indexes(prm_p(1:tr_p));
                
                test_number_per_cross_validation = min(tr_i, tr_p);
                
                interictal_test_index_subset = interictal_indexes(prm_i(end-test_number_per_cross_validation+1:end));
                preictal_test_index_subset = preictal_indexes(prm_p(end-test_number_per_cross_validation+1:end));
                
                train_index_subset = [interictal_train_index_subset ; preictal_train_index_subset];
                test_index_subset = [interictal_test_index_subset ; preictal_test_index_subset];
                
                r_train = r(train_index_subset, :);
                r_test = r(test_index_subset, :);
                r_all = r([train_index_subset ; test_index_subset], :);
                
                % train the classifier:
                % RBF SVM
                rbf_sigma = sigma(s) + p;
                box_constraint = f*length(interictal_train_index_subset)/length(preictal_train_index_subset);
                svmStruct = svmtrain(r_train, types(train_index_subset), 'kernel_function', 'rbf', 'rbf_sigma', rbf_sigma, 'boxconstraint', box_constraint);
                % svmStruct = svmtrain(r_train, types(train_index_subset), 'kernel_function', 'mlp');%, 'rbf_sigma', rbf_sigma, 'boxconstraint', length(interictal_train_index_subset)/length(preictal_train_index_subset));
                estimated_class = svmclassify(svmStruct, r_test);
                scores = double(estimated_class == preictal_label);
                [X, Y, T, AUC(k)] = perfcurve(types(test_index_subset), scores, preictal_label);
                estimated_class_test = svmclassify(svmStruct, r(test_indexes, :));
                PreictalPercentage(k) = 100*sum(double(estimated_class_test == preictal_label))/length(estimated_class_test);
            end
            str1 = ['Subject #', num2str(s), ', Sigma = ' num2str(rbf_sigma, '%5.2f') ', BCte = ' num2str(box_constraint, '%5.2f') ', AUC(min,mean,max)(SD) = (' num2str(min(AUC), '%5.2f') ',' num2str(mean(AUC), '%5.2f') ',' num2str(max(AUC), '%5.2f') ')(' num2str(std(AUC), '%5.2f') '), '];
            str2 = ['P%% = ' num2str(mean(PreictalPercentage), '%5.2f') '(' num2str(std(PreictalPercentage), '%5.2f') ')\n'];
            fprintf([str1 str2]);
            fprintf(fid, [str1 str2]);        
        end
    end
    fprintf('\n');
    % % %         fprintf(['I#:' num2str(length(interictal_indexes)) '\tP#:' num2str(length(preictal_indexes)) '\tT#:' num2str(length(test_indexes)) '\tI/P:' num2str(length(interictal_indexes)/length(preictal_indexes), '%5.2f') '\tpre: ' num2str(preictal_score, '%5.2f\t') '\tinter: ' num2str(interictal_score, '%5.2f\t'), '\ttest_ratio: ' num2str(test_ratio, '%5.2f\t') '\n']);
    % % %         for i = 1:length(test_indexes),
    % % %             fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(test_indexes(i)) > 0.5 );
    % % %             % % %         fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(test_indexes(i)));
    % % %             % % %         fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_vote_ratio(test_indexes(i)) - interictal_vote_ratio(test_indexes(i)) + 1.0);
    % % %         end
end
fclose(fid);