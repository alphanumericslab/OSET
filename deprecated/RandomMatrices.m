function [A, B, ang] = RandomMatrices(angle, dim1, dim2, flag)
% AdaptiveFilter has been deprecated. Use adaptive_filter instead.
warning('AdaptiveFilter has been deprecated. Use adaptive_filter instead.');
Kmax = 10000;
[A, B, ang] = random_mat_max_col_space_ang(angle, dim1, dim2, flag, Kmax);