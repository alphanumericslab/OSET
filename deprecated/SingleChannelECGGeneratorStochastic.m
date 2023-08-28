function [ecg, teta]= SingleChannelECGGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
% SingleChannelECGGeneratorStochastic has been deprecated. Use ecg_gen_stoochastic instead.
warning('SingleChannelECGGeneratorStochastic has been deprecated. Use ecg_gen_stoochastic instead.');
[ecg, teta]= ecg_gen_stoochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0);