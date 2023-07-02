% incorporate subject index in classification

clear
close all;
clc;

R1 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testFreqDomainEnergyIncreaseRate01.txt');
R2 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testICASubspaceAngles01.txt');
R3 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testWeightedHistogramAveraging02.txt');
R4 = importdata('/Users/rsameni/Desktop/Seizure/kaggle/SourcesWithResultsBackup/testProcessSeizureEEG_EnergyDistributions01.txt');

filenames = R1.textdata;
subject = R1.data(:, 1);
trial = R1.data(:, 2);
type = R1.data(:, 3);
r = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end)];
% r = R3.data(:, 4:end);

% discardcols = [];%[3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results
% r(:, discardcols) = [];

goodfeatures = [4 5 8 10 11 13 14 15 16 17 19 20 21 22 23 26 31 36 37 38];
r = r(:, goodfeatures);

sbj = 1;%(subject == 1);
column = 3;
bins = 75;

interictal_indexes = find(type == 1 & sbj); % interictal
preictal_indexes = find(type == 2 & sbj); % preictal
test_indexes = find(type == 3 & sbj); % test
train_indexes = [interictal_indexes ; preictal_indexes];
testfiles = filenames(test_indexes);

% per subject normalize columns
rr = zeros(size(r));
for s = 1:7,
    sb = (subject == s);
    rr(sb, :) = (r(sb, :) - ones(size(r(sb, :), 1), 1)*mean(r(sb, :), 1))./(ones(size(r(sb, :), 1), 1)*std(r(sb, :), [], 1));
end
% % % % % rr = r;

[n, r] = hist(rr(:, column), bins);
[n1, r1] = hist(rr(interictal_indexes, column), r);
[n2, r2] = hist(rr(preictal_indexes, column), r);
[n3, r3] = hist(rr(test_indexes, column), r);

N = sum(n);
N1 = sum(n1);
N2 = sum(n2);
N3 = sum(n3);


% % % rr = (rr - ones(size(rr,1),1)*mean(rr));
% % % rr = rr./(ones(size(rr,1),1)*std(rr));

% dropfeatures = 0; % 4
% Cx = cov(rr);
% [V, D] = eig(Cx);
% D = diag(D);
% [~, I] = sort(D, 1, 'descend');
% w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues
% 
% f = rr*w;
f = rr;

% b = TreeBagger(50, f(train_indexes, :), type(train_indexes) + 10*subject(train_indexes),'OOBPred','on','Method','classification');%,'prior', [0.5 0.5]);%, 'Cost',[0 1.0 ; 1.0 0]);
% [Y, scores, stdevs] = predict(b, f(test_indexes, :));
% figure
% plot(oobError(b));

% t = templateTree('minleaf',5);
t = ClassificationTree.template;
% ada = fitensemble(f(train_indexes, :), type(train_indexes), 'AdaBoostM1', 10000, t, 'prior', [0.5 5.0], 'type', 'classification');%%%%% ,  'method','Bag','type','classification');%,500,'tree');%,'RUSBoost',1000,'tree');%t,'LearnRate',0.1,'nprint',100);
ada = fitensemble(f(train_indexes, :), type(train_indexes) + 10*subject(train_indexes), 'AdaBoostM2', 100, t);%, 'prior', [0.5 5.0], 'type', 'classification');%%%%% ,  'method','Bag','type','classification');%,500,'tree');%,'RUSBoost',1000,'tree');%t,'LearnRate',0.1,'nprint',100);
[Y scores] = predict(ada, f(test_indexes, :));
figure
plot(resubLoss(ada,'mode','cumulative'));
grid

% 
% figure
% hold on
% h0 = bar(r, n/N, 'k');
% h1 = bar(r3, n3/N3, 'g');
% h2 = bar(r1, n1/N1, 'b');
% h3 = bar(r2, n2/N2, 'r');
% 
% p0 = findobj(h0,'Type','patch');
% p1 = findobj(h1,'Type','patch');
% p2 = findobj(h2,'Type','patch');
% p3 = findobj(h3,'Type','patch');
% % p4 = findobj(h4,'Type','patch');
% 
% set(p0, 'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
% set(p1, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.5);
% set(p2, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
% set(p3, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
% % set(p4, 'FaceColor', 'y', 'FaceAlpha', 0.5, 'EdgeColor', 'y', 'EdgeAlpha', 0.5);
% 
% legend('all', 'test', 'interictal', 'preictal');
% grid

% ldaClassFun= @(xtrain,ytrain,xtest)(classify(rr(test_indexes, :),rr(train_indexes, :), type(train_indexes)));
% ldaCVErr  = crossval('mcr',rr(train_indexes, :), type(train_indexes),'predfun', ...
%              ldaClassFun,'partition',cp);
         


% % % % % % % P = randperm(length(train_indexes)); % check to see of results depend on permutation of samples; no they didn't!
% % % 

% classifiertype = 'linear'; %'diaglinear'
% [class1, err1, POSTERIOR1, logp1, coeff1] = classify(f(test_indexes, :), f(train_indexes, :), type(train_indexes), classifiertype, [0.5 0.5]);


% % svmStruct = svmtrain(f(train_indexes, :), type(train_indexes), 'kernel_function', 'quadratic', 'kktviolationlevel', 0.15, 'method', 'QP');
% svmStruct = svmtrain(f(train_indexes, :), type(train_indexes), 'method', 'QP');
% 
% species = svmclassify(svmStruct, f(test_indexes, :));
% hold on;plot(5,2,'ro','MarkerSize',12);hold off

% % % % type = 'quadratic';
% % % % [class2, err2, POSTERIOR2, logp2, coeff2] = classify(f(test_indexes, :), f(train_indexes, :), type(train_indexes), type);
% % % 
% % % % % % % SVMStruct = svmtrain(f(train_indexes, :), type(train_indexes));
% % % % % % % Group = svmclassify(SVMStruct, f(test_indexes, :), 'Showplot', true);
% % % % % % %
% % % 

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

figure
plot(scores);
grid

figure
plot(Y);
grid

figure
plot(rr);
grid

% 
% fid = fopen('submission012.csv','w');
% fprintf(fid,'clip,preictal\n');
% for i = 1:length(species),
%     fprintf(fid,'%s,%10.8f\n', testfiles{i}, species(i)==2);
%     % % % %     fprintf(fid,'%s,%10.8f\n', testfiles{i}, min(max(scores(i, 2) - bias(2) + 0.5 - eps, 0), 1));
% end
% fclose(fid);
