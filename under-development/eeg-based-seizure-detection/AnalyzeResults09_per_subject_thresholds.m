clear
close all;
% clc;

R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');

filenames_all = R1.textdata;

subjects_all = R1.data(:, 1);

trials_all = R1.data(:, 2);

types_all = R1.data(:, 3);

r_all = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)])];

% r_all = [R5.data(:, [5:8 10:15 17:size(R5.data, 2)])];

 

% discardcols = [3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results

% r(:, discardcols) = [];

 

% goodfeatures = [4 5 8 10 11 13 14 15 16 17 19 20 21 22 23 26 31 36 37 38];

% r_all = r_all(:, goodfeatures);

 

fid = fopen('submission049.csv','w');

fprintf(fid,'clip,preictal\n');

 

% % % bins = [50 75 100 120 180 250 300 500 750];%75:6:81;%120;

% % % bins = 120;%50:1:550;%:13:250;%75:7:120;%75:6:81;%120;

% bins = 77;
% 
% % gmorder = 5;
% 
% dropfeatures = 5;

 

% % % % JUST FOR TEST

 

% per subject normalize columns

for s = 1:max(subjects_all),

    sb = (subjects_all == s);

    r_all(sb, :) = (r_all(sb, :) - ones(size(r_all(sb, :), 1), 1)*mean(r_all(sb, :), 1))./(ones(size(r_all(sb, :), 1), 1)*std(r_all(sb, :), [], 1));

    %         r_all(sb, :) = r_all(sb, :) - ones(size(r_all(sb, :), 1), 1)*mean(r_all(sb, :), 1);

end

 

% make the features uncorrelated
% 
% Cx = cov(r_all);
% 
% [V, D] = eig(Cx);
% 
% D = diag(D);
% 
% [~, I] = sort(D, 1, 'descend');
% 
% w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues
% 
% r_all = r_all*w;
% 
 

% sigma = [50 50 120 50 300 5 7.5];
sigma = [233 7 13 59 157 9 111];
bconst = [80 2.98 6.67 41.88 150 2.78 23.33];

for s = 1:max(subjects_all),

    subjects = (subjects_all == s);

    % % %     subjects = true(size(subjects_all)); % JUST FOR TEST

    

    filenames = filenames_all(subjects);

    trials = trials_all(subjects);

    types = types_all(subjects);

    rr = r_all(subjects, :);

    

    interictal_indexes = find(types == 1); % interictal

    preictal_indexes = find(types == 2); % preictal

    test_indexes = find(types == 3); % test

    train_indexes = [interictal_indexes ; preictal_indexes];%find(types == 1 | types == 2); % training set

    

    %testfiles = filenames(test_indexes);

    

    % % %     r_bins = cell(size(rr,2), length(bins));

    % % %     p = cell(size(rr,2), length(bins));

    % % %     p1 = cell(size(rr,2), length(bins));

    % % %     p2 = cell(size(rr,2), length(bins));

    % % %     p2_recovered = cell(size(rr,2), length(bins));

    % % %     p3 = cell(size(rr,2), length(bins));

    % % %     for feature = 1:size(rr,2),

    % % %         for k = 1:length(bins)

    % % %             [n, r_bins{feature, k}] = hist(rr(:, feature), bins(k));

    % % %

    % % %             [n1, r1] = hist(rr(interictal_indexes, feature), r_bins{feature, k});

    % % %             [n2, r2] = hist(rr(preictal_indexes, feature), r_bins{feature, k});

    % % %             [n3, r3] = hist(rr(test_indexes, feature), r_bins{feature, k});

    % % %             N = sum(n);

    % % %             N1 = sum(n1);

    % % %             N2 = sum(n2);

    % % %             N3 = sum(n3);

    % % %             p{feature, k} = n/N;

    % % %             p1{feature, k} = n1/N1;

    % % %             p2{feature, k} = n2/N2;

    % % %             p3{feature, k} = n3/N3;

    % % %

    % % %             p2_recovered{feature, k} = 2*(p3{feature, k} - 0.5*p1{feature, k}); % recover the preictal PDF from the test data assuming 50% interictal/preictal

    % % %             p2_recovered{feature, k}(p2_recovered{feature, k} < 0) = 0; % remove negative values

    % % %             p2_recovered{feature, k} = (length(test_indexes)/2*p2_recovered{feature, k} + length(preictal_indexes)*p2{feature, k})/(length(test_indexes)/2 + length(preictal_indexes));

    % % %             p2_recovered{feature, k} = p2_recovered{feature, k} / sum(p2_recovered{feature, k});

    % % %

    % % %             %             options = statset('MaxIter', 1000, 'TolFun', 1e-4);

    % % %             %             dstobj = gmdistribution.fit(rr(:, feature), gmorder, 'Options', options, 'Regularize', 1e-8);

    % % %             %             dstobj1 = gmdistribution.fit(rr(interictal_indexes, feature), gmorder, 'Options', options, 'Regularize', 1e-8);

    % % %             %             dstobj2 = gmdistribution.fit(rr(preictal_indexes, feature), gmorder, 'Options', options, 'Regularize', 1e-8);

    % % %             %             dstobj3 = gmdistribution.fit(rr(test_indexes, feature), gmorder, 'Options', options, 'Regularize', 1e-8);

    % % %             %             p{feature, k} = pdf(dstobj, r_bins{feature, k}');

    % % %             %             p1{feature, k} = pdf(dstobj1, r_bins{feature, k}');

    % % %             %             p2{feature, k} = pdf(dstobj2, r_bins{feature, k}');

    % % %             %             p3{feature, k} = pdf(dstobj3, r_bins{feature, k}');

    % % %             %             p{feature, k} = p{feature, k}/sum(p{feature, k});

    % % %             %             p1{feature, k} = p1{feature, k}/sum(p1{feature, k});

    % % %             %             p2{feature, k} = p2{feature, k}/sum(p2{feature, k});

    % % %             %             p3{feature, k} = p3{feature, k}/sum(p3{feature, k});

    % % %

    % % %             figure

    % % %             hold on

    % % %             h0 = bar(r_bins{feature, k}, p{feature, k}, 'k');

    % % %             h1 = bar(r_bins{feature, k}, p3{feature, k}, 'g');

    % % %             h2 = bar(r_bins{feature, k}, p1{feature, k}, 'b');

    % % %             h3 = bar(r_bins{feature, k}, p2_recovered{feature, k}, 'r');

    % % %             O0 = findobj(h0,'Type','patch');

    % % %             O1 = findobj(h1,'Type','patch');

    % % %             O2 = findobj(h2,'Type','patch');

    % % %             O3 = findobj(h3,'Type','patch');

    % % %             set(O0, 'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);

    % % %             set(O1, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.5);

    % % %             set(O2, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);

    % % %             set(O3, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);

    % % %             legend('all', 'test', 'interictal', 'preictal');

    % % %             % % % % %             legend('all', 'interictal');

    % % %             grid

    % % %         end

    % % %     end

    % % %

    % % %     interictal_vote = zeros(size(rr));

    % % %     preictal_vote = zeros(size(rr));

    % % %     for trial = 1:size(rr,1),

    % % %         for feature = 1:size(rr,2),

    % % %             for k = 1:length(bins)

    % % %                 I1 = find(r_bins{feature, k} > rr(trial, feature), 1, 'first');

    % % %                 I2 = max(I1 - 1, 1);

    % % %

    % % %                 % compare interictal with preictal probabilities

    % % %                 %                 delta = p1{feature, k}(I1) + p1{feature, k}(I2) - p2{feature, k}(I1) - p2{feature, k}(I2);

    % % %                 %                 if(~isempty(delta))

    % % %                 %                     if(delta > 0)

    % % %                 %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + 1;%delta;

    % % %                 %                     elseif(delta < 0)

    % % %                 %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + 1;%- delta;

    % % %                 %                     end

    % % %                 %                 end

    % % %

    % % %                 % compare interictal with preictal probabilities

    % % %                 delta = p1{feature, k}(I1) + p1{feature, k}(I2) - p2_recovered{feature, k}(I1) - p2_recovered{feature, k}(I2);

    % % %                 if(~isempty(delta))

    % % %                     if(delta > 0)

    % % %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + sum(abs(p1{feature, k} - p2_recovered{feature, k}))/2;%delta;

    % % %                     elseif(delta < 0)

    % % %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + sum(abs(p1{feature, k} - p2_recovered{feature, k}))/2;%- delta;

    % % %                     end

    % % %                 end

    % % %

    % % %                 % compare interictal with overall probabilities

    % % %                 %                 delta = p1{feature, k}(I1) + p1{feature, k}(I2) - p{feature, k}(I1) - p{feature, k}(I2);

    % % %                 %                 if(~isempty(delta))

    % % %                 %                     if(delta > 0)

    % % %                 %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + 1;%delta;%1;

    % % %                 %                     elseif(delta < 0)

    % % %                 %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + 1;%- delta;% 1;

    % % %                 %                     end

    % % %                 %                 end

    % % %

    % % %                 %                 % compare interictal with test and preictal probabilities

    % % %                 %                 delta_inter_vs_test = p1{feature, k}(I1) + p1{feature, k}(I2) - p3{feature, k}(I1) - p3{feature, k}(I2);

    % % %                 %                 delta_inter_vs_pre = p1{feature, k}(I1) + p1{feature, k}(I2) - p2{feature, k}(I1) - p2{feature, k}(I2);

    % % %                 %                 delta_test_vs_pre = p3{feature, k}(I1) + p3{feature, k}(I2) - p2{feature, k}(I1) - p2{feature, k}(I2);

    % % %                 %                 if(~isempty(delta_inter_vs_test) && ~isempty(delta_inter_vs_pre))

    % % %                 %                     if(delta_inter_vs_test > 0 && delta_inter_vs_pre > 0)

    % % %                 %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + 1;%sum(abs(p1{feature, k} - p3{feature, k}))/2;

    % % %                 %                     end

    % % %                 %                 end

    % % %                 %                 if(~isempty(delta_inter_vs_test) && ~isempty(delta_inter_vs_pre))

    % % %                 %                     if(delta_inter_vs_test < 0 && delta_inter_vs_pre < 0)

    % % %                 %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + 1;%sum(abs(p1{feature, k} - p3{feature, k}))/2;

    % % %                 %                     end

    % % %                 %                 end

    % % %

    % % %                 %                 % compare interictal with test probabilities (not a good idea; because all test data are classified as preictal)

    % % %                 %                 delta_inter_vs_test = p1{feature, k}(I1) + p1{feature, k}(I2) - p3{feature, k}(I1) - p3{feature, k}(I2);

    % % %                 %                 % % % guard = 0.001;

    % % %                 %                 if(~isempty(delta_inter_vs_test))

    % % %                 %                     if(delta_inter_vs_test > 0)

    % % %                 %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + sum(abs(p1{feature, k} - p3{feature, k}))/2;

    % % %                 %                     elseif(delta_inter_vs_test < 0)

    % % %                 %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + sum(abs(p1{feature, k} - p3{feature, k}))/2;

    % % %                 %                     end

    % % %                 %                 end

    % % %

    % % %             end

    % % %         end

    % % %     end

    % % %

    % % %     preictal_vote_ratio = sum(preictal_vote, 2)./(sum(preictal_vote, 2) + sum(interictal_vote, 2));

    % % %     interictal_vote_ratio = sum(interictal_vote, 2)./(sum(preictal_vote, 2) + sum(interictal_vote, 2));

    

    % test other classifiers:

    % % %     b = TreeBagger(100, rr(train_indexes, :), types(train_indexes),'OOBPred','on','Method','classification','MinLeaf', 1, 'prior', [0.01 10.0]);%, 'NVarToSample', 10);%, 'Cost',[0 1.0 ; 1.0 0]);

    % % %     figure

    % % %     plot(oobError(b));

    % % %     [Y, pp, stdevs] = predict(b, rr);

    

    % % % % % t = templateTree('minleaf',5);

    % % %     rusTree = fitensemble(rr(train_indexes, :), types(train_indexes),'AdaBoostM1',100,'tree');%t,'LearnRate',0.1,'nprint',100);

    % % %     [Y, pp] = predict(rusTree, rr);

    

    % % %     cvens = crossval(rusTree);%, 'kfold', 100);

    % % %     L = 100*kfoldLoss(cvens)

    

    % % %     interictal_vote_ratio = pp(:, 1);

    % % %     preictal_vote_ratio = pp(:, 2);

    

    % % %     train_indexes_balanced = [interictal_indexes ; repmat(preictal_indexes, ceil(length(interictal_indexes)/length(preictal_indexes)), 1)];

    % % %     svmStruct = svmtrain(rr(train_indexes_balanced, :), types(train_indexes_balanced), 'kernel_function', 'rbf');

    %     svmStruct = svmtrain(rr(train_indexes, :), types(train_indexes), 'kernel_function', 'rbf', 'boxconstraint', [length(interictal_indexes)*ones(length(interictal_indexes),1) ; length(preictal_indexes)*ones(length(preictal_indexes),1)]);

    % % % std(rr)

    % % % length(interictal_indexes)/length(preictal_indexes)

    

    % Single SVM

    %     svmStruct = svmtrain(rr(train_indexes, :), types(train_indexes), 'kernel_function', 'rbf', 'rbf_sigma', sigma(s), 'boxconstraint', length(interictal_indexes)/length(preictal_indexes));% std(std(rr))%, 'boxconstraint', length(interictal_indexes)/length(preictal_indexes));%[length(interictal_indexes)*ones(length(interictal_indexes),1) ; length(preictal_indexes)*ones(length(preictal_indexes),1)]);%, 'method', 'QP');%, 'boxconstraint', .1*ones(length(train_indexes),1));

    %     Y = svmclassify(svmStruct, rr);

    %     I_inter = find(Y == 1);

    %     I_pre = find(Y == 2);

    %     interictal_vote_ratio = zeros(length(rr), 1);   interictal_vote_ratio(I_inter) = 1;

    %     preictal_vote_ratio = zeros(length(rr), 1);     preictal_vote_ratio(I_pre) = 1;

    

    % Voted SVM

    interictal_vote_ratio = zeros(length(rr), 1);

    preictal_vote_ratio = zeros(length(rr), 1);

%     sigma2 = [50 50 120 50 300 5 7.5];

% % %     sigma2 = [50 27 30 50 300 5 6];

 

% % %     for d = -0:.6:0,

        %         for e = -.2:.2:.2,

%         svmStruct = svmtrain(rr(train_indexes, :), types(train_indexes), 'kernel_function', 'rbf', 'rbf_sigma', sigma2(s)+d, 'boxconstraint', length(interictal_indexes)/length(preictal_indexes));% std(std(rr))%, 'boxconstraint', length(interictal_indexes)/length(preictal_indexes));%[length(interictal_indexes)*ones(length(interictal_indexes),1) ; length(preictal_indexes)*ones(length(preictal_indexes),1)]);%, 'method', 'QP');%, 'boxconstraint', .1*ones(length(train_indexes),1));
        svmStruct = svmtrain(rr(train_indexes, :), types(train_indexes), 'kernel_function', 'rbf', 'rbf_sigma', sigma(s), 'boxconstraint', bconst(s));

        Y = svmclassify(svmStruct, rr);

        I_inter = find(Y == 1);

        I_pre = find(Y == 2);

        interictal_vote_ratio(I_inter) = interictal_vote_ratio(I_inter) + 1;

        preictal_vote_ratio(I_pre) = preictal_vote_ratio(I_pre) + 1;

        %         end

% % %     end

    

    score = (sign( (preictal_vote_ratio - interictal_vote_ratio) .* (2*types-3)) + 1)/2.0;

    preictal_prob = (preictal_vote_ratio - interictal_vote_ratio + 1.0)/2.0;
%     preictal_prob = double(Y == 2);

    

    preictal_score = 100*sum(score(preictal_indexes))/length(preictal_indexes);

    interictal_score = 100*sum(score(interictal_indexes))/length(interictal_indexes);

    test_ratio = 100*sum(preictal_vote_ratio(test_indexes) > interictal_vote_ratio(test_indexes) ) / length(test_indexes);

    

    

    fprintf(['I#:' num2str(length(interictal_indexes)) '\tP#:' num2str(length(preictal_indexes)) '\tT#:' num2str(length(test_indexes)) '\tI/P:' num2str(length(interictal_indexes)/length(preictal_indexes), '%5.2f') '\tpre: ' num2str(preictal_score, '%5.2f\t') '\tinter: ' num2str(interictal_score, '%5.2f\t'), '\ttest_ratio: ' num2str(test_ratio, '%5.2f\t') '\n']);

    

    figure

    subplot(121);

    plot(2*types(train_indexes)-3, 'r', 'linewidth', 2.5);

    hold on;

    stem(preictal_vote_ratio(train_indexes) - interictal_vote_ratio(train_indexes));

    grid

    title('train');

    a = axis;

    

    subplot(122);

    stem(preictal_vote_ratio(test_indexes) - interictal_vote_ratio(test_indexes), 'r');

    grid

    title('test');

    aa = axis;

    aa(3) = a(3);

    aa(4) = a(4);

    axis(aa);

    

    for i = 1:length(test_indexes),

        fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(test_indexes(i)) > 0.5 );

        % % %         fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(test_indexes(i)));

        % % %         fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_vote_ratio(test_indexes(i)) - interictal_vote_ratio(test_indexes(i)) + 1.0);

    end

end

% % fclose(fid);