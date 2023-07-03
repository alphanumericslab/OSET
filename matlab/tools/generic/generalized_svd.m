function [U, S, V] = generalized_svd(A, X, Y)
%
% [U, S, V] = generalized_svd(A, X, Y)
% Generalized Singular Value Decomposition
%
% input:
% A: recrangular matrix (m x n)
% X: square positive definite matrix (m x m); row-wise constraints on A
% Y: square positive definite matrix (n x n); column-wise constraints on A
%
% outputs:
% U, S, and V, such that A = U*S*V' and U'*X*U = Im and V'*Y*V = In
%
% Reference:
%   Abdi, Hervé. "Singular value decomposition (SVD) and generalized singular
%       value decomposition." Encyclopedia of measurement and statistics (2007): 
%       907-912. URL: https://www.utdallas.edu/~herve/Abdi-SVD2007-pretty.pdf
%
%   Revision History:
%       2019: First release
%       2023: Renamed from deprecated version GeneralizedSVD()
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

x = sqrtm(X);
y = sqrtm(Y);

B = x * A * y;
[Q, S, R] = svd(B);

U = x\Q;
V = y\R;
