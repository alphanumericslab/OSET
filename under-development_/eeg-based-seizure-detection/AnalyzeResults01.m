% my first submission resulted in a score: 0.57797

clear
close all;

R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');

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
dropfeatures = 4; % 4

% normalize columns
rr = (r - ones(size(r, 1), 1)*mean(r, 1))./(ones(size(r, 1), 1)*std(r, [], 1));

Cx = cov(rr);

[V, D] = eig(Cx);
D = diag(D);
[Y,I] = sort(D, 1, 'descend');
w = V(:, I(1:end - dropfeatures)); % remove the last few zero eigenvalues

f = rr*w;

interictal_indexes = find(mode == 1); % interictal
preictal_indexes = find(mode == 2); % preictal
test_indexes = find(mode == 3); % test
train_indexes = [interictal_indexes ; preictal_indexes];

% P = randperm(length(train_indexes)); % check to see of results depend on permutation of samples; no they didn't!

testfiles = filenames(test_indexes);

type = 'linear'; %'diaglinear'
[class1, err1, POSTERIOR1, logp1, coeff1] = classify(f(test_indexes, :), f(train_indexes, :), mode(train_indexes), type);
type = 'quadratic';
[class2, err2, POSTERIOR2, logp2, coeff2] = classify(f(test_indexes, :), f(train_indexes, :), mode(train_indexes), type);

figure
plot(class1-class2, 'ro')
grid

figure
hold on
plot(class1-1,'ko');
plot(POSTERIOR1);
grid

fid = fopen('submission002.csv','w');
fprintf(fid,'clip,preictal\n');
for i = 1:length(POSTERIOR2),
    fprintf(fid,'%s,%8.6f\n', testfiles{i}, POSTERIOR2(i, 2));
end
fclose(fid);

% stem(D(I));
% grid
% r = data;
% % bins = 0:.01:1;
% % bins = -10:.5:50;
% bins = 30;%0:1:100;
% column = 4;
%
% subject = 1;%(r(:, 1) == 5);
% mode1 = (r(:, 3) == 1); % interictal
% mode2 = (r(:, 3) == 2); % preictal
% mode3 = (r(:, 3) == 3); % test
%
% [n1, r1] = hist(r(mode1 & subject, column), bins);
% [n2, r2] = hist(r(mode2 & subject, column), bins);
% [n3, r3] = hist(r(mode3 & subject, column), bins);
%
% N1 = sum(n1);
% N2 = sum(n2);
% N3 = sum(n3);
%
% figure
% hold on
% bar(r1, n1/N1);
% bar(r2, n2/N2, 'r');
% % bar(r3, n3/N3, 'g');
%
% % h = findobj(gca,'Type','patch');
% % set(h(1), 'FaceColor', 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'EdgeAlpha', 0.9);
% % set(h(2), 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.9);
% % set(h(3), 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.9);
% % legend('interictal', 'preictal', 'test');
% grid
