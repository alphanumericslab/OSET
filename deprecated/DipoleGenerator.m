function [DIP, teta] = DipoleGenerator(N, fs, f, alphai, bi, tetai, teta0)
% DipoleGenerator has been deprecated. Use ecg_dipole_gen_state_space instead.
warning('DipoleGenerator has been deprecated. Use ecg_dipole_gen_state_space instead.');
[DIP, teta] = ecg_dipole_gen_state_space(N, fs, f, alphai, bi, tetai, teta0);