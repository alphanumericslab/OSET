function [DIP, teta] = DipoleGenerator(N, fs, f, alphai, bi, tetai, teta0)
% DipoleGenerator has been deprecated. Use vcg_gen_state_space instead.
warning('DipoleGenerator has been deprecated. Use vcg_gen_state_space instead.');
[DIP, teta] = vcg_gen_state_space(N, fs, f, alphai, bi, tetai, teta0);