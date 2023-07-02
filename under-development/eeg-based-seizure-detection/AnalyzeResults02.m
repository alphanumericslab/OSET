clear
close all;

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
mode = R1.data(:, 3);
r = [R1.data(:, 4:end) R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end)];
% r = R3.data(:, 4:end);
discardcols = [];%[3 24 25 28 29 30 31 34 35 36]; % submission004.csv; didn't improve the results

r(:, discardcols) = [];

dropfeatures = 4; % 4
interictal_indexes = find(mode == 1); % interictal
preictal_indexes = find(mode == 2); % preictal
test_indexes = find(mode == 3); % test
train_indexes = [interictal_indexes ; preictal_indexes];
testfiles = filenames(test_indexes);

column = 1;
bins = 100;
[n1, r1] = hist(r(interictal_indexes, column), bins);
[n2, r2] = hist(r(preictal_indexes, column), bins);
[n3, r3] = hist(r(test_indexes, column), bins);

N1 = sum(n1);
N2 = sum(n2);
N3 = sum(n3);

figure
hold on
% bar(r3, n3/N3, 'g');
bar(r1, n1/N1, 'b');
bar(r2, n2/N2, 'r');
grid

% normalize columns
rr = (r - ones(size(r, 1), 1)*mean(r, 1))./(ones(size(r, 1), 1)*std(r, [], 1));

Cx = cov(rr);

[V, D] = eig(Cx);
D = diag(D);
[~, I] = sort(D, 1, 'descend');
w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues

f = rr*w;

% P = randperm(length(train_indexes)); % check to see of results depend on permutation of samples; no they didn't!

type = 'linear'; %'diaglinear'
[class1, err1, POSTERIOR1, logp1, coeff1] = classify(f(test_indexes, :), f(train_indexes, :), mode(train_indexes), type);
type = 'quadratic';
[class2, err2, POSTERIOR2, logp2, coeff2] = classify(f(test_indexes, :), f(train_indexes, :), mode(train_indexes), type);

% SVMStruct = svmtrain(f(train_indexes, :), mode(train_indexes));
% Group = svmclassify(SVMStruct, f(test_indexes, :), 'Showplot', true);

b = TreeBagger(50, f(train_indexes, :), mode(train_indexes),'OOBPred','on','Method','classification','prior', [0.5 0.5]);%, 'Cost',[0 1.0 ; 1.0 0]);
[Y, scores, stdevs] = predict(b, f(test_indexes, :));

%t = templateTree('minleaf',5);
% rusTree = fitensemble(f(train_indexes, :), mode(train_indexes),'RUSBoost',1000,'tree');%t,'LearnRate',0.1,'nprint',100);
% [Y, scores, stdevs] = predict(rusTree, f(test_indexes, :));


sum(str2num(cell2mat(Y))-1)/length(test_indexes)

figure
stem(str2num(cell2mat(Y))-1.5);
grid

figure
plot(scores);
grid

figure
plot(stdevs);
grid

figure
plot(oobError(b))
xlabel('number of grown trees');
ylabel('out-of-bag classification error');


figure
plot(class1-class2, 'ro')
grid

figure
hold on
plot(class1-1,'ko');
plot(POSTERIOR1);
grid

fid = fopen('submission006.csv','w');
fprintf(fid,'clip,preictal\n');
for i = 1:length(POSTERIOR1),
%    fprintf(fid,'%s,%10.8f\n', testfiles{i}, POSTERIOR1(i, 2));
    fprintf(fid,'%s,%10.8f\n', testfiles{i}, min(scores(i, 2) + 0.439,1));
end
fclose(fid);

% % stem(D(I));
% % grid
% % r = data;
% % % bins = 0:.01:1;
% % % bins = -10:.5:50;
% % bins = 30;%0:1:100;
% % column = 4;
% %
% % subject = 1;%(r(:, 1) == 5);
% % mode1 = (r(:, 3) == 1); % interictal
% % mode2 = (r(:, 3) == 2); % preictal
% % mode3 = (r(:, 3) == 3); % test
% %
% % [n1, r1] = hist(r(mode1 & subject, column), bins);
% % [n2, r2] = hist(r(mode2 & subject, column), bins);
% % [n3, r3] = hist(r(mode3 & subject, column), bins);
% %
% % N1 = sum(n1);
% % N2 = sum(n2);
% % N3 = sum(n3);
% %
% % figure
% % hold on
% % bar(r1, n1/N1);
% % bar(r2, n2/N2, 'r');
% % % bar(r3, n3/N3, 'g');
% %
% % % h = findobj(gca,'Type','patch');
% % % set(h(1), 'FaceColor', 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'EdgeAlpha', 0.9);
% % % set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.9);
% % % set(h(3), 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.9);
% % % legend('interictal', 'preictal', 'test');
% % grid
