function [CosTheta, U, V] = subspace_principal_angles(A, B)
% subspace_principal_angles calculates principal angles and vectors between two column subspaces.
%
% Syntax: [CosTheta, U, V] = subspace_principal_angles(A, B)
% 
% Inputs:
%   A: Matrix representing the first column subspace (size: m x n1, where m is the number of rows and n1 is the number of columns).
%   B: Matrix representing the second column subspace (size: m x n2, where m is the number of rows and n2 is the number of columns).
%
% Outputs:
%   CosTheta: Principal angles between the two subspaces (size: min(n1, n2) x 1).
%   U: Matrix containing the principal vectors associated with the first subspace (size: m x min(n1, n2)).
%   V: Matrix containing the principal vectors associated with the second subspace (size: m x min(n1, n2)).
%
% Reference:
%   Golub, G.H. & Van Loan, C.F. (1996). Matrix Computations, Section 12.4.3, Page 603.
%
%   Revision History:
%       2022: First release
%       2023: Renamed from deprecated version SubspaceAngles()
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Check input sizes
if size(A, 1) ~= size(B, 1)
    error('Column space size mismatch; A and B should have the same number of rows');
end

if size(A, 2) < size(B, 2)
    error('Number of columns of B should not exceed the number of columns of A');
end

% Compute QR factorization of A and B
[Qa, ~] = qr(A, 0);
[Qb, ~] = qr(B, 0);

% Calculate the matrix C and perform singular value decomposition on it
C = Qa' * Qb;

[Y, S, Z] = svd(C);

% Compute principal angles
S = min(diag(S), 1);
CosTheta = max(S, -1);

% Compute principal vectors
U = Qa * Y;
V = Qb * Z;
