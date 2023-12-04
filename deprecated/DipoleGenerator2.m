function [DIP, teta] = DipoleGenerator2(N, fs, f, alphai, bi, tetai, teta0)
% DipoleGenerator2 has been deprecated. Use vcg_gen_direct_sum instead.
warning('DipoleGenerator2 has been deprecated. Use vcg_gen_direct_sum instead.');
[DIP, teta] = vcg_gen_direct_sum(N, fs, f, alphai, bi, tetai, teta0);