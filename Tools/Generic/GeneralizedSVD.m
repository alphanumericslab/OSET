function [U, S, V] = GeneralizedSVD(A, X, Y)
%
% [U, S, V] = GeneralizedSVD(A, X, Y)
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
% Abdi, Hervé. "Singular value decomposition (SVD) and generalized singular
% value decomposition." Encyclopedia of measurement and statistics (2007): 
% 907-912.
% URL: https://www.utdallas.edu/~herve/Abdi-SVD2007-pretty.pdf
%
% Open Source Electrophysiological Toolbox, version 3.14, June 2019
% Released under the GNU General Public License
% Copyright (C) 2019  Reza Sameni
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

x = sqrtm(X);
y = sqrtm(Y);

B = x * A * y;
[Q, S, R] = svd(B);

U = x\Q;
V = y\R;
