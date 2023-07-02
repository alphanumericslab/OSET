clear
close all;
clc;

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

fid = fopen('submission041.csv','w');
fprintf(fid,'clip,preictal\n');

% bins = [50 75 100 120 180 250 300 500 750];%75:6:81;%120;
bins = 120;%:13:250;%75:7:120;%75:6:81;%120;
gmorder = 5;

% per subject normalize columns
% JUST FOR TEST
for s = 1:max(subjects_all),
    sb = (subjects_all == s);
    r_all(sb, :) = (r_all(sb, :) - ones(size(r_all(sb, :), 1), 1)*mean(r_all(sb, :), 1))./(ones(size(r_all(sb, :), 1), 1)*std(r_all(sb, :), [], 1));
    %     r_all(sb, :) = r_all(sb, :) - ones(size(r_all(sb, :), 1), 1)*mean(r_all(sb, :), 1);
end


for s = 1:max(subjects_all),
    subjects = (subjects_all == s);
    % % % % % %     subjects = true(size(subjects_all)); % JUST FOR TEST
    
    filenames = filenames_all(subjects);
    trials = trials_all(subjects);
    types = types_all(subjects);
    rr = r_all(subjects, :);
    
    interictal_indexes = find(types == 1); % interictal
    preictal_indexes = find(types == 2); % preictal
    test_indexes = find(types == 3); % test
    train_indexes = [interictal_indexes ; preictal_indexes];%find(types == 1 | types == 2); % training set
    
    %testfiles = filenames(test_indexes);
    
    r_bins = cell(size(rr,2), length(bins));
    p = cell(size(rr,2), length(bins));
    p1 = cell(size(rr,2), length(bins));
    p2 = cell(size(rr,2), length(bins));
    p3 = cell(size(rr,2), length(bins));
    for feature = 1:size(rr,2),
        for k = 1:length(bins)
            [n, r_bins{feature, k}] = hist(rr(:, feature), bins(k));
            
            [n1, r1] = hist(rr(interictal_indexes, feature), r_bins{feature, k});
            [n2, r2] = hist(rr(preictal_indexes, feature), r_bins{feature, k});
            [n3, r3] = hist(rr(test_indexes, feature), r_bins{feature, k});
            N = sum(n);
            N1 = sum(n1);
            N2 = sum(n2);
            N3 = sum(n3);
            p{feature, k} = n/N;
            p1{feature, k} = n1/N1;
            p2{feature, k} = n2/N2;
            p3{feature, k} = n3/N3;
            
            %             options = statset('MaxIter', 1000, 'TolFun', 1e-4);
            %             dstobj = gmdistribution.fit(rr(:, feature), gmorder, 'Options', options, 'Regularize', 1e-8);
            %             dstobj1 = gmdistribution.fit(rr(interictal_indexes, feature), gmorder, 'Options', options, 'Regularize', 1e-8);
            %             dstobj2 = gmdistribution.fit(rr(preictal_indexes, feature), gmorder, 'Options', options, 'Regularize', 1e-8);
            %             dstobj3 = gmdistribution.fit(rr(test_indexes, feature), gmorder, 'Options', options, 'Regularize', 1e-8);
            %             p{feature, k} = pdf(dstobj, r_bins{feature, k}');
            %             p1{feature, k} = pdf(dstobj1, r_bins{feature, k}');
            %             p2{feature, k} = pdf(dstobj2, r_bins{feature, k}');
            %             p3{feature, k} = pdf(dstobj3, r_bins{feature, k}');
            %             p{feature, k} = p{feature, k}/sum(p{feature, k});
            %             p1{feature, k} = p1{feature, k}/sum(p1{feature, k});
            %             p2{feature, k} = p2{feature, k}/sum(p2{feature, k});
            %             p3{feature, k} = p3{feature, k}/sum(p3{feature, k});
            
            %             figure
            %             hold on
            %             h0 = bar(r_bins{feature, k}, p{feature, k}, 'k');
            %             %         h1 = bar(r_bins{feature, k}, p3{feature, k}, 'g');
            %             h2 = bar(r_bins{feature, k}, p1{feature, k}, 'b');
            %             %         h3 = bar(r_bins{feature, k}, p2{feature, k}, 'r');
            %             O0 = findobj(h0,'Type','patch');
            %             %         O1 = findobj(h1,'Type','patch');
            %             O2 = findobj(h2,'Type','patch');
            %             %         O3 = findobj(h3,'Type','patch');
            %             %             set(O0, 'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
            %             %             %         set(O1, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.5);
            %             %             set(O2, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
            %             %             %         set(O3, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
            %             %             %         legend('all', 'test', 'interictal', 'preictal');
            %             legend('all', 'interictal');
            %             grid
        end
    end
    
    interictal_vote = zeros(size(rr));
    preictal_vote = zeros(size(rr));
    for trial = 1:size(rr,1),
        for feature = 1:size(rr,2),
            for k = 1:length(bins)
                I1 = find(r_bins{feature, k} > rr(trial, feature), 1, 'first');
                I2 = max(I1 - 1, 1);
                
                % compare interictal with preictal probabilities
                %                 delta = p1{feature, k}(I1) + p1{feature, k}(I2) - p2{feature, k}(I1) - p2{feature, k}(I2);
                %                 if(~isempty(delta))
                %                     if(delta > 0)
                %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + 1;%delta;
                %                     elseif(delta < 0)
                %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + 1;%- delta;
                %                     end
                %                 end
                
                % compare interictal with overall probabilities
                %                 delta = p1{feature, k}(I1) + p1{feature, k}(I2) - p{feature, k}(I1) - p{feature, k}(I2);
                %                 if(~isempty(delta))
                %                     if(delta > 0)
                %                         interictal_vote(trial, feature) = interictal_vote(trial, feature) + 1;%delta;%1;
                %                     elseif(delta < 0)
                %                         preictal_vote(trial, feature) = preictal_vote(trial, feature) + 1;%- delta;% 1;
                %                     end
                %                 end
                
                % compare interictal with test probabilities
                delta_inter_vs_test = p1{feature, k}(I1) + p1{feature, k}(I2) - p3{feature, k}(I1) - p3{feature, k}(I2);
                delta_inter_vs_pre = p1{feature, k}(I1) + p1{feature, k}(I2) - p2{feature, k}(I1) - p2{feature, k}(I2);
                delta_test_vs_pre = p3{feature, k}(I1) + p3{feature, k}(I2) - p2{feature, k}(I1) - p2{feature, k}(I2);
                if(~isempty(delta_inter_vs_test) && ~isempty(delta_inter_vs_pre))
                    if(delta_inter_vs_test > 0 && delta_inter_vs_pre > 0)
                        interictal_vote(trial, feature) = interictal_vote(trial, feature) + 1;%sum(abs(p1{feature, k} - p3{feature, k}))/2;
                    end
                end
                if(~isempty(delta_inter_vs_test) && ~isempty(delta_inter_vs_pre))
                    if(delta_inter_vs_test < 0 && delta_inter_vs_pre < 0)
                        preictal_vote(trial, feature) = preictal_vote(trial, feature) + 1;%sum(abs(p1{feature, k} - p3{feature, k}))/2;
                    end
                end
            end
        end
    end
    
    preictal_vote_ratio = sum(preictal_vote, 2)./(sum(preictal_vote, 2) + sum(interictal_vote, 2));
    interictal_vote_ratio = sum(interictal_vote, 2)./(sum(preictal_vote, 2) + sum(interictal_vote, 2));

    % test other classifiers:
b = TreeBagger(150, rr(train_indexes, :), types(train_indexes),'OOBPred','on','Method','classification','prior', [0.1 0.9]);%, 'Cost',[0 1.0 ; 1.0 0]);
[Y, pp, stdevs] = predict(b, rr);
interictal_vote_ratio = pp(:, 1);
preictal_vote_ratio = pp(:, 2);
% end test other classifiers

    score = (sign( (preictal_vote_ratio - interictal_vote_ratio) .* (2*types-3)) + 1)/2.0;
    preictal_prob = (preictal_vote_ratio - interictal_vote_ratio + 1.0)/2.0;
    


% % % % % % % t = templateTree('minleaf',5);
% % % % % % % rusTree = fitensemble(f(train_indexes, :), type(train_indexes),'RUSBoost',1000,'tree');%t,'LearnRate',0.1,'nprint',100);
% % % % % % % [Y, scores, stdevs] = predict(rusTree, f(test_indexes, :));
% % % % % % %
    
    preictal_score = 100*sum(score(preictal_indexes))/length(preictal_indexes);
    interictal_score = 100*sum(score(interictal_indexes))/length(interictal_indexes);
    test_ratio = 100*sum(preictal_vote_ratio(test_indexes) > interictal_vote_ratio(test_indexes) ) / length(test_indexes);
    
    
    disp(['pre: ' num2str(preictal_score) ', inter: ' num2str(interictal_score), ', test_ratio: ' num2str(test_ratio)]);
    
    figure
    subplot(121);
    plot(2*types(train_indexes)-3, 'r', 'linewidth', 2.5);
    hold on;
    stem(preictal_vote_ratio(train_indexes) - interictal_vote_ratio(train_indexes));
    grid
    title('train');
    
    subplot(122);
    stem(preictal_vote_ratio(test_indexes) - interictal_vote_ratio(test_indexes), 'r');
    grid
    title('test');
    
    for i = 1:length(test_indexes),
%         fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(test_indexes(i)) > 0.5 );
        fprintf(fid,'%s,%10.8f\n', filenames{test_indexes(i)}, preictal_prob(test_indexes(i)) );
    end
end
% % fclose(fid);
