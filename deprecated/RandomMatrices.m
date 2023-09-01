function [A, B, ang] = RandomMatrices(angle, dim1, flag)
% RandomMatrices has been deprecated. Use random_mat_max_col_space_ang instead.
warning('RandomMatrices has been deprecated. Use random_mat_max_col_space_ang instead.');
Kmax = 10000;
[A, B, ang] = random_mat_max_col_space_ang(angle, dim1, 3, flag, Kmax);