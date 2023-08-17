function [dipole, teta]= DipoleGeneratorStochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0)
% DipoleGeneratorStochastic has been deprecated. Use ecg_dipole_gen_stochastic instead.
warning('DipoleGeneratorStochastic has been deprecated. Use ecg_dipole_gen_stochastic instead.');
[dipole, teta]= ecg_dipole_gen_stochastic(N, fs, f, f_deviations, alphai, delta_alphai, bi, delta_bi, tetai, delta_tetai, teta0);