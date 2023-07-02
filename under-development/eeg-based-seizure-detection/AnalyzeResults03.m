clear
close all;
clc;

pth = 'E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\';
R1 = importdata([pth 'testFreqDomainEnergyIncreaseRate01.txt']);
R2 = importdata([pth 'testICASubspaceAngles01.txt']);
R3 = importdata([pth 'testWeightedHistogramAveraging02.txt']);
R4 = importdata([pth 'testProcessSeizureEEG_EnergyDistributions01.txt']);
R5 = importdata([pth 'testSpectralAnalysis01.txt']);

% r1 = R1.data;
% r2 = R2.data;
% r3 = R3.data;
% r4 = R4.data;

filenames = R1.textdata;
subject = R1.data(:, 1);
trial = R1.data(:, 2);
type = R1.data(:, 3);
r = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end)];
r = R5.data(:, 4:end);
discardcols = [];%[3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results

r(:, discardcols) = [];

sbj = 1;%(subject == 6);
column = 1;
bins = 100;

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

rr = r;

[n1, r1] = hist(rr(interictal_indexes, column), bins);
[n2, r2] = hist(rr(preictal_indexes, column), bins);
[n3, r3] = hist(rr(test_indexes, column), bins);

N1 = sum(n1);
N2 = sum(n2);
N3 = sum(n3);

figure
hold on
h1 = bar(r3, n3/N3, 'g');
h2 = bar(r1, n1/N1, 'b');
h3 = bar(r2, n2/N2, 'r');

p1 = findobj(h1,'Type','patch');
p2 = findobj(h2,'Type','patch');
p3 = findobj(h3,'Type','patch');
set(p1, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.5);
set(p2, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'EdgeAlpha', 0.5);
set(p3, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'EdgeAlpha', 0.5);
legend('test', 'interictal', 'preictal');
grid

% % % dropfeatures = 4; % 4
% % % 
% % % Cx = cov(rr);
% % % 
% % % [V, D] = eig(Cx);
% % % D = diag(D);
% % % [~, I] = sort(D, 1, 'descend');
% % % w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues
% % % 
% % % f = rr*w;
% % % 
% % % % % % % P = randperm(length(train_indexes)); % check to see of results depend on permutation of samples; no they didn't!
% % % 
% % % classifiertype = 'linear'; %'diaglinear'
% % % [class1, err1, POSTERIOR1, logp1, coeff1] = classify(f(test_indexes, :), f(train_indexes, :), type(train_indexes), classifiertype);
% % % 
% % % % type = 'quadratic';
% % % % [class2, err2, POSTERIOR2, logp2, coeff2] = classify(f(test_indexes, :), f(train_indexes, :), type(train_indexes), type);
% % % 
% % % % % % % SVMStruct = svmtrain(f(train_indexes, :), type(train_indexes));
% % % % % % % Group = svmclassify(SVMStruct, f(test_indexes, :), 'Showplot', true);
% % % % % % %
% % % 
% % % b = TreeBagger(50, f(train_indexes, :), type(train_indexes),'OOBPred','on','Method','classification','prior', [0.5 0.5]);%, 'Cost',[0 1.0 ; 1.0 0]);
% % % [Y, scores, stdevs] = predict(b, f(test_indexes, :));
% % % 
% % % % % % % t = templateTree('minleaf',5);
% % % % % % % rusTree = fitensemble(f(train_indexes, :), type(train_indexes),'RUSBoost',1000,'tree');%t,'LearnRate',0.1,'nprint',100);
% % % % % % % [Y, scores, stdevs] = predict(rusTree, f(test_indexes, :));
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
% % % hold on
% % % plot(class1-1,'ko');
% % % plot(POSTERIOR1);
% % % grid
% % % 
% % % fid = fopen('submission009.csv','w');
% % % fprintf(fid,'clip,preictal\n');
% % % for i = 1:length(scores),
% % %     fprintf(fid,'%s,%10.8f\n', testfiles{i}, POSTERIOR1(i, 2));
% % %     % % % %     fprintf(fid,'%s,%10.8f\n', testfiles{i}, min(max(scores(i, 2) - bias(2) + 0.5 - eps, 0), 1));
% % % end
% % % fclose(fid);
% % % 
% % % % % % % stem(D(I));
% % % % % % % grid
% % % % % % % r = data;
% % % % % % % % bins = 0:.01:1;
% % % % % % % % bins = -10:.5:50;
% % % % % % % bins = 30;%0:1:100;
% % % % % % % column = 4;
% % % % % % %
% % % % % % % subject = 1;%(r(:, 1) == 5);
% % % % % % % type1 = (r(:, 3) == 1); % interictal
% % % % % % % type2 = (r(:, 3) == 2); % preictal
% % % % % % % type3 = (r(:, 3) == 3); % test
% % % % % % %
% % % % % % % [n1, r1] = hist(r(type1 & subject, column), bins);
% % % % % % % [n2, r2] = hist(r(type2 & subject, column), bins);
% % % % % % % [n3, r3] = hist(r(type3 & subject, column), bins);
% % % % % % %
% % % % % % % N1 = sum(n1);
% % % % % % % N2 = sum(n2);
% % % % % % % N3 = sum(n3);
% % % % % % %
% % % % % % % figure
% % % % % % % hold on
% % % % % % % bar(r1, n1/N1);
% % % % % % % bar(r2, n2/N2, 'r');
% % % % % % % % bar(r3, n3/N3, 'g');
% % % % % % %
