function [theta_min, theta_max, theta_mean, theta_med, ang] = ColumnSubspaceAngle(A, B)
% Deprecated: ColumnSubspaceAngle is deprecated. Use column_subspace_angle instead.
warning('Deprecated: ColumnSubspaceAngle is deprecated. Use column_subspace_angle instead.');
[ang, theta_min, theta_max, theta_mean, theta_med] = column_subspace_angle(A, B);