clear
close all;
clc;

R1 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testFreqDomainEnergyIncreaseRate01.txt');
R2 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testICASubspaceAngles01.txt');
R3 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testWeightedHistogramAveraging02.txt');
R4 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testProcessSeizureEEG_EnergyDistributions01.txt');

% r1 = R1.data;
% r2 = R2.data;
% r3 = R3.data;
% r4 = R4.data;

filenames = R1.textdata;
subject = R1.data(:, 1);
trial = R1.data(:, 2);
type = R1.data(:, 3);
r = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end)];

% discardcols = [3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results
% r(:, discardcols) = [];

goodfeatures = [4 5 8 10 11 13 14 15 16 17 19 20 21 22 23 26 31 36 37 38];
r = r(:, goodfeatures);

sbj = (subject == 4);% | (subject == 2);
bins = 75;%75:10:120;
gmorder = 10;

interictal_indexes = find(type == 1 & sbj); % interictal
preictal_indexes = find(type == 2 & sbj); % preictal
test_indexes = find(type == 3 & sbj); % test
train_indexes = [interictal_indexes ; preictal_indexes];
testfiles = filenames(test_indexes);

% per subject normalize columns
% rr = zeros(size(r));
% for s = 1:7,
%     sb = (subject == s);
%     rr(sb, :) = (r(sb, :) - ones(size(r(sb, :), 1), 1)*mean(r(sb, :), 1))./(ones(size(r(sb, :), 1), 1)*std(r(sb, :), [], 1));
%     %     rr(sb, :) = (r(sb, :) - ones(size(r(sb, :), 1), 1)*mean(r(sb, :), 1));
% end
rr = r;

r_bins = cell(size(rr,2), length(bins));
p = cell(size(rr,2), length(bins));
p1 = cell(size(rr,2), length(bins));
p2 = cell(size(rr,2), length(bins));
p3 = cell(size(rr,2), length(bins));
for column = 1:size(rr,2),
    for k = 1:length(bins)
        [n, r_bins{column, k}] = hist(rr(:, column), bins(k));
        
        [n1, r1] = hist(rr(interictal_indexes, column), r_bins{column, k});
        [n2, r2] = hist(rr(preictal_indexes, column), r_bins{column, k});
        [n3, r3] = hist(rr(test_indexes, column), r_bins{column, k});
        N = sum(n);
        N1 = sum(n1);
        N2 = sum(n2);
        N3 = sum(n3);
        p{column, k} = n/N;
        p1{column, k} = n1/N1;
        p2{column, k} = n2/N2;
        p3{column, k} = n3/N3;
        
        %         options = statset('MaxIter', 5000);
        %         dstobj = gmdistribution.fit(rr(:, column), gmorder, 'Options', options);
        %         dstobj1 = gmdistribution.fit(rr(interictal_indexes, column), gmorder, 'Options', options);
        %         dstobj2 = gmdistribution.fit(rr(preictal_indexes, column), gmorder, 'Options', options);
        %         dstobj3 = gmdistribution.fit(rr(test_indexes, column), gmorder, 'Options', options);
        %         p{column, k} = pdf(dstobj, r_bins{column, k}');
        %         p1{column, k} = pdf(dstobj1, r_bins{column, k}');
        %         p2{column, k} = pdf(dstobj2, r_bins{column, k}');
        %         p3{column, k} = pdf(dstobj3, r_bins{column, k}');
        %         p{column, k} = p{column, k}/sum(p{column, k});
        %         p1{column, k} = p1{column, k}/sum(p1{column, k});
        %         p2{column, k} = p2{column, k}/sum(p2{column, k});
        %         p3{column, k} = p3{column, k}/sum(p3{column, k});
        
        figure
        hold on
        h0 = bar(r_bins{column, k}, p{column, k}, 'k');
        %         h1 = bar(r_bins{column, k}, p3{column, k}, 'g');
        h2 = bar(r_bins{column, k}, p1{column, k}, 'b');
        %         h3 = bar(r_bins{column, k}, p2{column, k}, 'r');
        O0 = findobj(h0,'Type','patch');
        %         O1 = findobj(h1,'Type','patch');
        O2 = findobj(h2,'Type','patch');
        %         O3 = findobj(h3,'Type','patch');
        set(O0, 'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
        %         set(O1, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.5);
        set(O2, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
        %         set(O3, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
        %         legend('all', 'test', 'interictal', 'preictal');
        legend('all', 'interictal');
        
        grid
    end
end

% prob_preictal = zeros(size(rr,1));
interictal_vote = zeros(size(rr));
preictal_vote = zeros(size(rr));
for trial = 1:size(rr,1),
    for column = 1:size(rr,2),
        for k = 1:length(bins)
            I1 = find(r_bins{column, k} > rr(trial, column), 1, 'first');
            I2 = max(I1 - 1, 1);
            
            % compare interictal with preictal probabilities
            delta = p1{column, k}(I1) + p1{column, k}(I2) - p2{column, k}(I1) - p2{column, k}(I2);
            if(~isempty(delta))
                if(delta > 0)
                    % if( double(p1{column, k}(I1) + p1{column, k}(I2)) > double(p2{column, k}(I1) + p2{column, k}(I2) ))
                    interictal_vote(trial, column) = interictal_vote(trial, column) + 1;%delta;
                else
                    preictal_vote(trial, column) = preictal_vote(trial, column) + 1;%- delta;
                end
            end
            
            % compare interictal with overall probabilities
            %             if(p1{column, k}(I1) + p1{column, k}(I2)> p{column, k}(I1) + p{column, k}(I2))
            %                 interictal_vote(trial, column) = interictal_vote(trial, column) + 1;
            %             else
            %                 preictal_vote(trial, column) = preictal_vote(trial, column) + 1;
            %             end
            
            % compare interictal with test probabilities
            %             if(p1{column, k}(I1) + p1{column, k}(I2)> p3{column, k}(I1) + p3{column, k}(I2))
            %                 interictal_vote(trial, column) = interictal_vote(trial, column) + 1;
            %             else
            %                 preictal_vote(trial, column) = preictal_vote(trial, column) + 1;
            %             end
            
        end
    end
end

%     prob1 = sum(interictal_vote) / (sum(interictal_vote) + sum(preictal_vote) );
%     prob2 = sum(preictal_vote) / (sum(interictal_vote) + sum(preictal_vote) );
%     if(prob1 > prob2)
%         disp([num2str(trial) ': interictal, p = ' num2str(prob1)]);
%     else
%         disp([num2str(trial) ': preictal,   p = ' num2str(prob2)]);
%     end
%     prob_preictal(trial) = prob2;

prob_preictal = preictal_vote - interictal_vote;

% tp1 = 2*type(interictal_indexes)-3;
% e1 = sum(abs(prob_preictal(interictal_indexes,:) - tp1*ones(1,size(prob_preictal(interictal_indexes,:), 2))), 1)/length(interictal_indexes)
% tp2 = 2*type(preictal_indexes)-3;
% e2 = sum(abs(prob_preictal(preictal_indexes,:) - tp2*ones(1,size(prob_preictal(preictal_indexes,:), 2))), 1)/length(preictal_indexes)


figure
plot(2*type(train_indexes)-3, 'r', 'linewidth', 2.5);
hold on;
stem(mean(prob_preictal(train_indexes, :), 2));
grid

% % % % ldaClassFun= @(xtrain,ytrain,xtest)(classify(rr(test_indexes, :),rr(train_indexes, :), type(train_indexes)));
% % % % ldaCVErr  = crossval('mcr',rr(train_indexes, :), type(train_indexes),'predfun', ...
% % % %              ldaClassFun,'partition',cp);
% % %
% % %
% % % rr = (rr - ones(size(rr,1),1)*mean(rr));
% % % rr = rr./(ones(size(rr,1),1)*std(rr));
% % %
% % % dropfeatures = 4; % 4
% % % Cx = cov(rr);
% % % [V, D] = eig(Cx);
% % % D = diag(D);
% % % [~, I] = sort(D, 1, 'descend');
% % % w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues
% % %
% % % f = rr*w;
% % %
% % % % % % % % % % P = randperm(length(train_indexes)); % check to see of results depend on permutation of samples; no they didn't!
% % % % % %
% % %
% % % classifiertype = 'linear'; %'diaglinear'
% % % [class1, err1, POSTERIOR1, logp1, coeff1] = classify(f(test_indexes, :), f(train_indexes, :), type(train_indexes), classifiertype, [0.5 0.5]);
% % %
% % %
% % % % % svmStruct = svmtrain(f(train_indexes, :), type(train_indexes), 'kernel_function', 'quadratic', 'kktviolationlevel', 0.15, 'method', 'QP');
% % % % svmStruct = svmtrain(f(train_indexes, :), type(train_indexes), 'method', 'QP');
% % % %
% % % % species = svmclassify(svmStruct, f(test_indexes, :));
% % % % hold on;plot(5,2,'ro','MarkerSize',12);hold off
% % %
% % % % % % % type = 'quadratic';
% % % % % % % [class2, err2, POSTERIOR2, logp2, coeff2] = classify(f(test_indexes, :), f(train_indexes, :), type(train_indexes), type);
% % % % % %
% % % % % % % % % % SVMStruct = svmtrain(f(train_indexes, :), type(train_indexes));
% % % % % % % % % % Group = svmclassify(SVMStruct, f(test_indexes, :), 'Showplot', true);
% % % % % % % % % %
% % % % % %
% % % % % % b = TreeBagger(50, f(train_indexes, :), type(train_indexes),'OOBPred','on','Method','classification','prior', [0.5 0.5]);%, 'Cost',[0 1.0 ; 1.0 0]);
% % % % % % [Y, scores, stdevs] = predict(b, f(test_indexes, :));
% % % % % %
% % % % % % % % % % t = templateTree('minleaf',5);
% % % % % % % % % % rusTree = fitensemble(f(train_indexes, :), type(train_indexes),'RUSBoost',1000,'tree');%t,'LearnRate',0.1,'nprint',100);
% % % % % % % % % % [Y, scores, stdevs] = predict(rusTree, f(test_indexes, :));
% % % % % % %
% % %
% % % sum(str2num(cell2mat(Y))-1)/length(test_indexes)
% % %
% % % bias = median(scores, 1)
% % %
% % % figure
% % % stem(str2num(cell2mat(Y))-1.5);
% % % grid
% % %
% % % figure
% % % plot(scores);
% % % grid
% % %
% % % figure
% % % plot(stdevs);
% % % grid
% % %
% % % figure
% % % plot(oobError(b))
% % % xlabel('number of grown trees');
% % % ylabel('out-of-bag classification error');
% % %
% % %
% % % % % % figure
% % % % % % plot(class1-class2, 'ro')
% % % % % % grid
% % %

% % % figure
% % % % hold on
% % % % plot(class1-1,'ko');
% % % % plot(POSTERIOR1);
% % % % plot(min(max(prob_preictal - median(prob_preictal) + 0.5, 0), 1));
% % % plot(prob_preictal);
% % % grid

% fid = fopen('submission017.csv','w');
% fprintf(fid,'clip,preictal\n');
% for i = 1:length(prob_preictal),
%     fprintf(fid,'%s,%10.8f\n', testfiles{i}, double(prob_preictal(i)));
%     %     fprintf(fid,'%s,%10.8f\n', testfiles{i}, min(max(prob_preictal(i)-median(prob_preictal) + 0.5,0),1));
%     % % % %     fprintf(fid,'%s,%10.8f\n', testfiles{i}, min(max(scores(i, 2) - bias(2) + 0.5 - eps, 0), 1));
% end
% fclose(fid);
