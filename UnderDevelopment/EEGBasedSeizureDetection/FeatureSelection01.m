clear
close all;
clc

per_subject = 1;
outputfileame = 'Submission50.csv';
outfid = fopen(outputfileame, 'w');
fprintf(outfid,'clip,preictal\n');

% load data
R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');
% R6 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures3.txt');

filenames_all = R1.textdata;
subjects_all = R1.data(:, 1);
trials_all = R1.data(:, 2);
types_all = R1.data(:, 3);
interictal_label = 1;
preictal_label = 2;
test_label = 3;
r_raw = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)])];% R6.data(:, 4:end)];

fun = @(XT,yT,Xt,yt)...
    (sum(~strcmp(yt, classify(Xt,XT,yT,'quadratic'))));

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
    
    c = cvpartition(types(train_indexes), 'k', 5);
    %     opts = statset('UseParallel','always');
    
    [featureset, history] = sequentialfs(fun, r(train_indexes, :), types(train_indexes), 'cv', c);%, 'options', opts);
    %     [featureset, history] = sequentialfs(@FSCost, r(train_indexes, :), types(train_indexes), 'cv', c);%,'options',opts);
    
    %     svmStruct = svmtrain(r(train_indexes, logical(featureset)), types(train_indexes), 'kernel_function', 'rbf', 'rbf_sigma', 100, 'boxconstraint', 5);
    %     estimated_class = svmclassify(svmStruct, r(:, logical(featureset)));
    
    [estimated_class error] = classify(r(:, logical(featureset)), r(train_indexes, logical(featureset)), types(train_indexes), 'quadratic');
    
    error
    
    figure
    subplot(211);
    stem(featureset);
    grid
    title(num2str(s));
    
    subplot(223);
    stem(estimated_class(train_indexes) - 1.5);
    grid
    
    subplot(224);
    stem(estimated_class(test_indexes) - 1.5, 'r');
    grid
    
    for i = 1:length(test_indexes),
        fprintf(outfid,'%s,%10.8f\n', filenames{test_indexes(i)}, double(estimated_class(test_indexes(i)) == preictal_label));
    end
    
end
fclose(outfid);
