% load data
% % % R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
% % % R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
% % % R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
% % % R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
% % % R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');
% % % R6 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures3.txt');

% R = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeatures_phase.txt');
% AnalyzeResults14_SpectralFeaturesFunction(R.data, R.textdata, 'testProcessSeizureEEGPhaseAmplitudeFeatures_phase.csv', 'testProcessSeizureEEGPhaseAmplitudeFeatures_phase_unbiased.csv', 1, 850);
% disp('#1 Done!');

% R = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_amp.txt');
% AnalyzeResults14_SpectralFeaturesFunction(R.data, R.textdata, 'testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_amp.csv', 'testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_amp_unbiased.csv', 1, 850);
% disp('#2 Done!');
% 
% R = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeatures_amp.txt');
% AnalyzeResults14_SpectralFeaturesFunction(R.data, R.textdata, 'testProcessSeizureEEGPhaseAmplitudeFeatures_amp.csv', 'testProcessSeizureEEGPhaseAmplitudeFeatures_amp_unbiased.csv', 1, 850);
% disp('#3 Done!');
% 
% R = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_phase.txt');
% AnalyzeResults14_SpectralFeaturesFunction(R.data, R.textdata, 'testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_phase.csv', 'testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_phase_unbiased.csv', 1, 850);
% disp('#4 Done!');

% R = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp_wavelet.txt');
% AnalyzeResults14_SpectralFeaturesFunction(R.data, R.textdata, 'testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp_wavelet.csv', 'testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp_wavelet_unbiased.csv', 1, 850);
% disp('#5 Done!');

% R = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp.txt');
% AnalyzeResults14_SpectralFeaturesFunction(R.data, R.textdata, 'testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp.csv', 'testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp_unbiased.csv', 1, 850);
% disp('#6 Done!');

% R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
% R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
% R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
% R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
% R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');
% R6 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures3.txt');
% % R7 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising_amp_unbiased.txt');
% % R8 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGPhaseAmplitudeFeaturesXCorrBased_amp_wavelet.txt');
% 
% Rdata = [R1.data R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)]) R6.data(:, 4:end)];% R7.data(:, 4:end) R8.data(:, 4:end)];
% Rtextdata = R1.textdata;
% disp('Data loaded!');
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures1000.csv', 'AllFeatures1000_unbiased.csv', 1, 1000);
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures1000NoPerSubject.csv', 'AllFeatures1000NoPerSubject_unbiased.csv', 0, 1000);
% 
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures2000.csv', 'AllFeatures2000_unbiased.csv', 1, 2000);
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures2000NoPerSubject.csv', 'AllFeatures2000NoPerSubject_unbiased.csv', 0, 2000);
% 
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures3000.csv', 'AllFeatures3000_unbiased.csv', 1, 3000);
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures3000NoPerSubject.csv', 'AllFeatures3000NoPerSubject_unbiased.csv', 0, 3000);
% 
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures5000.csv', 'AllFeatures5000_unbiased.csv', 1, 5000);
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'AllFeatures5000NoPerSubject.csv', 'AllFeatures5000NoPerSubject_unbiased.csv', 0, 5000);
% disp('Done!');
% 
% R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGSpectralWaveletFeatures4.txt');
% R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures4.txt');
% Rdata = [R1.data R2.data(:, 4:end)];
% Rtextdata = R1.textdata;
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'testProcessSeizureEEGSpectralAndWaveletFeatures4.csv', 'testProcessSeizureEEGSpectralAndWaveletFeatures4_unbiased.csv', 1, 1000);
% disp('Done!');

% R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
% R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
% R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
% R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
% R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');
% R6 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures3.txt');
% R7 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGRWAofRawSpectra.txt');
% 
% Rdata = [R1.data R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)]) R6.data(:, 4:end) R7.data(:, 4:end)];% R8.data(:, 4:end)];
% Rtextdata = R1.textdata;
% disp('Data loaded!');
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'testProcessSeizureEEGRWAofRawSpectra1000.csv', 'testProcessSeizureEEGRWAofRawSpectra1000_unbiased.csv', 1, 1000);
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'testProcessSeizureEEGRWAofRawSpectra850.csv', 'testProcessSeizureEEGRWAofRawSpectra850_unbiased.csv', 1, 850);
% AnalyzeResults14_SpectralFeaturesFunction(Rdata, Rtextdata, 'testProcessSeizureEEGRWAofRawSpectra3000.csv', 'testProcessSeizureEEGRWAofRawSpectra3000_unbiased.csv', 1, 3000);
% disp('Done!');

% % % R7 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGRWAofRawSpectra.txt');
% % % Rdata = R7.data;
% % % Rtextdata = R7.textdata;

% R1 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testFreqDomainEnergyIncreaseRate01.txt');
% R2 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testICASubspaceAngles01.txt');
% R3 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testWeightedHistogramAveraging02.txt');
% R4 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEG_EnergyDistributions01.txt');
% R5 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGICAOnNormalizedHistogram.txt');
R6 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGAllSpectralFeatures3.txt');
R7 = importdata('E:\Sameni\Projects\Seizure\SourceCodes\Probabilistic\SourcesWithResultsBackup\testProcessSeizureEEGRWAofRawSpectra.txt');

Rdata = [R6.data R7.data(:, 4:end)];%[R1.data R2.data(:, 4:end) R3.data(:, 4:end) R4.data(:, 4:end) R5.data(:, [5:8 10:15 17:size(R5.data, 2)]) R6.data(:, 4:end) R7.data(:, 4:end)];% R8.data(:, 4:end)];
Rtextdata = R6.textdata;
disp('Data loaded!');

% AnalyzeResults15_14Class(Rdata, Rtextdata, 'test1000.csv', 'test_unbiased.csv', 1000);
% AnalyzeResults16_SpectralFeaturesFunction(Rdata, Rtextdata, 'testLast1000NoPriors.csv', 'testLast1000NoPriors_unbiased.csv', 1, 1000);
AnalyzeResults17_SpectralFeaturesFunction(Rdata, Rtextdata, 'test_post_contest_linear_classifier_lastfew_features.csv', 'test_post_contest_linear_classifier_lastfew_features_unbiased.csv', 1);
disp('Done!');