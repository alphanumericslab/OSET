function [CosTheta, U, V] = SubspaceAngles(A, B)
% Calculate principal angles and vectors between two column subspaces
%
% Reza Sameni (reza.sameni@gmail.com)
% Feb 2022
%
% Ref: G.H. Golub & C.F. Van Loan (1996), Matrix Computations, Section 12.4.3, Page 603

if size(A, 1) ~= size(B, 1)
    error('Column space size mis-match; A and B should have the same number of rows');
end

if size(A, 2) < size(B, 2)
    error('Number of columns of B should not exceed the number of columns of A');
end

[Qa, ~] = qr(A, 0);
[Qb, ~] = qr(B, 0);
C = Qa' * Qb;
[Y, S, Z] = svd(C);
S = min(diag(S), 1);
CosTheta = max(S, -1);
U = Qa * Y;
V = Qb * Z;
