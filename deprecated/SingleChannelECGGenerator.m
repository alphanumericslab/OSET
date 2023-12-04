function [ECG, teta] = SingleChannelECGGenerator(teta, teta0, alphai, bi, tetai)
% SingleChannelECGGenerator has been deprecated. Use ecg_gen_gmm instead.
warning('SingleChannelECGGenerator has been deprecated. Use ecg_gen_gmm instead.');
[ECG, teta] = ecg_gen_gmm(teta, teta0, alphai, bi, tetai);