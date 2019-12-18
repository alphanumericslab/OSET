close all;
clear;
clc;

% testFreqDomainEnergyIncreaseRate01
% 
% testICASubspaceAngles01
% 
% testWeightedHistogramAveraging02
% 
% testProcessSeizureEEG_EnergyDistributions01

% testProcessSeizureEEGAllSpectralFeatures
% 
% AnalyzeResults11_WithCrossValidationResults

% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction350.csv', 'AnalyzeResults13_SpectralFeaturesFunction350_unbiased.csv', 1, 350);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction450.csv', 'AnalyzeResults13_SpectralFeaturesFunction450_unbiased.csv', 1, 450);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction550.csv', 'AnalyzeResults13_SpectralFeaturesFunction550_unbiased.csv', 1, 550);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction650.csv', 'AnalyzeResults13_SpectralFeaturesFunction650_unbiased.csv', 1, 650);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction750.csv', 'AnalyzeResults13_SpectralFeaturesFunction750_unbiased.csv', 1, 750);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction825.csv', 'AnalyzeResults13_SpectralFeaturesFunction825_unbiased.csv', 1, 825);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction850.csv', 'AnalyzeResults13_SpectralFeaturesFunction850_unbiased.csv', 1, 850);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction875.csv', 'AnalyzeResults13_SpectralFeaturesFunction875_unbiased.csv', 1, 875);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction900.csv', 'AnalyzeResults13_SpectralFeaturesFunction900_unbiased.csv', 1, 900);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction925.csv', 'AnalyzeResults13_SpectralFeaturesFunction925_unbiased.csv', 1, 925);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction950.csv', 'AnalyzeResults13_SpectralFeaturesFunction950_unbiased.csv', 1, 950);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction975.csv', 'AnalyzeResults13_SpectralFeaturesFunction975_unbiased.csv', 1, 975);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1050.csv', 'AnalyzeResults13_SpectralFeaturesFunction1050_unbiased.csv', 1, 1050);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1150.csv', 'AnalyzeResults13_SpectralFeaturesFunction1150_unbiased.csv', 1, 1150);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1250.csv', 'AnalyzeResults13_SpectralFeaturesFunction1250_unbiased.csv', 1, 1250);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1350.csv', 'AnalyzeResults13_SpectralFeaturesFunction1350_unbiased.csv', 1, 1350);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1450.csv', 'AnalyzeResults13_SpectralFeaturesFunction1450_unbiased.csv', 1, 1450);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1550.csv', 'AnalyzeResults13_SpectralFeaturesFunction1550_unbiased.csv', 1, 1550);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1650.csv', 'AnalyzeResults13_SpectralFeaturesFunction1650_unbiased.csv', 1, 1650);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1750.csv', 'AnalyzeResults13_SpectralFeaturesFunction1750_unbiased.csv', 1, 1750);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1850.csv', 'AnalyzeResults13_SpectralFeaturesFunction1850_unbiased.csv', 1, 1850);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction1950.csv', 'AnalyzeResults13_SpectralFeaturesFunction1950_unbiased.csv', 1, 1950);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction2050.csv', 'AnalyzeResults13_SpectralFeaturesFunction2050_unbiased.csv', 1, 2050);
% AnalyzeResults13_SpectralFeaturesFunction('AnalyzeResults13_SpectralFeaturesFunction2150.csv', 'AnalyzeResults13_SpectralFeaturesFunction2150_unbiased.csv', 1, 2150);

testProcessSeizureEEGPhaseAmplitudeFeatures_AfterDenoising

testProcessSeizureEEGPhaseAmplitudeFeatures

