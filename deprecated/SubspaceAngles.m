function [CosTheta, U, V] = SubspaceAngles(A, B)
% SubspaceAngles has been deprecated. Use subspace_principal_angles instead.
warning('SubspaceAngles has been deprecated. Use subspace_principal_angles instead.');
[CosTheta, U, V] = subspace_principal_angles(A, B);