clear
close all;
clc

% parameters
logfileame = 'SVMOptimizedParamters03.txt';
outputfileame = 'SVMOptimizedParamters03Submission.csv';
% % % bins = [50 75 100 120 180 250 300 500 750];%75:6:81;%120;
% % % bins = 120;%50:1:550;%:13:250;%75:7:120;%75:6:81;%120;
% % % bins = 77;
% gmorder = 5;
% dropfeatures = 5;
sigma = 150*ones(1, 7);
% sigma = [181.5 27.0 179 66 111 26 19];
per_subject = 1;
cross_validation_number = 5;
train_ratio_per_cross_validation = 0.5; % the ratio of the interictal/preictal data used for test for each cross-validation iteration
minAUCThreshold = 0.53;
% meanAUCThreshold = 0.6;
% maxAUCThreshold = 0.75;
minPreictalPercentage = 40;
maxPreictalPercentage = 60;

% fix the random seed for randperm
rng(0);

% load data
R1 = importdata('testProcessSeizureEEGAllSpectralFeatures3.txt');
% % % R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
% % % R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
% % % R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
% % % R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
% % % R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');

filenames_all = R1.textdata;
subjects_all = R1.data(:, 1);
trials_all = R1.data(:, 2);
types_all = R1.data(:, 3);
interictal_label = 1;
preictal_label = 2;
test_label = 3;
dropfeatures = 0;
% % % r_raw = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)])];
r_raw = R1.data(:, 4:end);

% discard/keep some features
% discardcols = [3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results
% r_raw(:, discardcols) = [];
% goodfeatures = [4 5 8 10 11 13 14 15 16 17 19 20 21 22 23 26 31 36 37 38];
% r_raw = r_raw(:, goodfeatures);

% log file
logfid = fopen(logfileame, 'w');
outfid = fopen(outputfileame, 'w');
fprintf(outfid,'clip,preictal\n');

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
% % %
% % % figure
% % % stem(D(I)/D(I(1)));
% % % grid

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
    train_indexes = [interictal_indexes ; preictal_indexes];
    test_indexes = find(types == test_label); % test
    %testfiles = filenames(test_indexes);
    
    preictal_votes = zeros(length(test_indexes), 1);
    interictal_votes = zeros(length(test_indexes), 1);
    number_of_votes = 0;
    % cross validate over labeled data
    for f = [1./(10:-1:1) (2:10)],%[1/10 1/9 1/8 1/5 1/4 1/3 1/2 1 2 3 4 5 10],
        for p = -149:149,
            AUC = zeros(cross_validation_number, 1);
            tr_i = floor(train_ratio_per_cross_validation*length(interictal_indexes));
            tr_p = floor(train_ratio_per_cross_validation*length(preictal_indexes));
            rbf_sigma = sigma(s) + p;
            box_constraint = f*tr_i/tr_p;
            for k = 1:cross_validation_number,
                prm_i = randperm(length(interictal_indexes));
                interictal_train_index_subset = interictal_indexes(prm_i(1:tr_i));
                
                prm_p = randperm(length(preictal_indexes));
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
                svmStruct = svmtrain(r_train, types(train_index_subset), 'kernel_function', 'rbf', 'rbf_sigma', rbf_sigma, 'boxconstraint', box_constraint);
                estimated_class = svmclassify(svmStruct, r_test);
                scores = double(estimated_class == preictal_label);
                [X, Y, T, AUC(k)] = perfcurve(types(test_index_subset), scores, preictal_label);
            end
            minAUC = min(AUC);
            meanAUC = mean(AUC);
            maxAUC = max(AUC);
            stdAUC = std(AUC);
            % if(minAUC > minAUCThreshold || meanAUC > meanAUCThreshold || maxAUC > maxAUCThreshold)
            if(minAUC > minAUCThreshold)
                svmStruct = svmtrain(r(train_indexes, :), types(train_indexes), 'kernel_function', 'rbf', 'rbf_sigma', rbf_sigma, 'boxconstraint', box_constraint);
                estimated_class_test = svmclassify(svmStruct, r(test_indexes, :));
                PreictalPercentage = 100*sum(double(estimated_class_test == preictal_label))/length(estimated_class_test);
                if(PreictalPercentage > minPreictalPercentage && PreictalPercentage < maxPreictalPercentage)
                    I_i = (estimated_class_test == interictal_label);
                    I_p = (estimated_class_test == preictal_label);
                    interictal_votes(I_i) = interictal_votes(I_i) + (meanAUC - 0.5);
                    preictal_votes(I_p) = preictal_votes(I_p) + (meanAUC - 0.5);
                    number_of_votes = number_of_votes + 1;
                    str1 = ['Subject #', num2str(s), ', Vote # ' num2str(number_of_votes) ', Sigma = ' num2str(rbf_sigma, '%5.2f') ', BCte = ' num2str(box_constraint, '%5.2f') ', AUC(min,mean,max)(SD) = (' num2str(minAUC, '%5.2f') ',' num2str(meanAUC, '%5.2f') ',' num2str(maxAUC, '%5.2f') ')(' num2str(stdAUC, '%5.2f') '), '];
                    str2 = ['P%% = ' num2str(PreictalPercentage, '%5.2f') '\n'];
                    fprintf([str1 str2]);
                    fprintf(logfid, [str1 str2]);
                end
            end
        end
    end
    if(number_of_votes > 0)
        preictal_prob = (preictal_votes - interictal_votes)/number_of_votes + 0.5;
    else
        preictal_prob = rand(length(test_indexes), 1);  % zeros(length(test_indexes), 1) + 0.5;
    end
    for i = 1:length(test_indexes),
        fprintf(outfid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(i));
    end
end
fclose(logfid);
fclose(outfid);