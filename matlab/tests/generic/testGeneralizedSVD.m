% test function for Generalized Singular Value Decomposition
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

clear;
clc;
close all;

m = 3; % row dimensions
n = 5; % column dimensions

% the matrix used in SVD; replace with the desired matrix
A = randn(m, n)

% the PD matrix used as the row-wise constraints on A; replace with the
% desired constraint matrix
x = randn(m);
X = x*x'
X_eigs = eig(X)

% the PD matrix used as the column-wise constraints on A; replace with the
% desired constraint matrix
y = randn(n);
Y = y*y'
Y_eigs = eig(Y)

[U, S, V] = GeneralizedSVD(A, X, Y)

reconstruction_error = A - U*S*V'
row_constraint = U'*X*U
column_constraint = V'*Y*V