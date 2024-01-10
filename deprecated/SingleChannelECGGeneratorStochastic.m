function [ecg, teta]= SingleChannelECGGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
% SingleChannelECGGeneratorStochastic has been deprecated. Use ecg_gen_stochastic instead.
warning('SingleChannelECGGeneratorStochastic has been deprecated. Use ecg_gen_stochastic instead.');
[ecg, teta]= ecg_gen_stochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0);