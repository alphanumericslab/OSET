function [theta_min, theta_max, theta_mean, theta_med, ang] = ColumnSubspaceAngle(A, B)
% ColumnSubspaceAngle has been deprecated. Use column_subspace_angle instead.
warning('ColumnSubspaceAngle has been deprecated. Use column_subspace_angle instead.');
[ang, theta_min, theta_max, theta_mean, theta_med] = column_subspace_angle(A, B);