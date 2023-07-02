% An implementation of Jacobi's algorithm for symmetric matrix
% diagonalization: The classical Jacobi algorithm.
% Ref: Golub, Gene H., and Charles F. Van Loan. Matrix computations.
% JHU Press, 1996, Sec 8.4
%
% Reza Sameni 3/1/2019
% email: reza.sameni@gmail.com
%
function testJacobiDiagonalization02
clc;
clear;
close all;

% A = [6 5 1; 5 22 3 ; 1 3 7]; % counter example in which the algorithm does not converge
% % A = [9 5 1; 5 8 3 ; 1 3 7];
% N = 3;

A = [6 5 -1 4; 5 22 3 2; -1 3 37 5 ; 4 2 5 10]; % counter example in which the algorithm does not converge
N = 4;

ITR = 25;
P = A
e = zeros(1, ITR);
U = eye(N);
for k = 1 : ITR
    B = triu(abs(P), 1);
    [i, j] = find(B == max(B(:)), 1);
    [c, s] = SymSchur2(P(i, i), P(j, j), P(i, j));
    Uk = eye(N);
    Uk(i, i) = c;
    Uk(i, j) = s;
    Uk(j, i) = -s;
    Uk(j, j) = c;
    
    P = Uk' * P * Uk
    e(k) = sum(sum((P - diag(diag(P))).^2));
    U = U * Uk
end

[V,D] = eig(A)

figure;
plot(1:ITR, 10*log10(e));
grid;
title('The log diagonalization error profile per iteration');

function [c s] = SymSchur2(App, Aqq, Apq)
if(Apq ~= 0)
    tau = (Aqq - App)/(2*Apq);
    if(tau >= 0)
        t = 1/(tau + sqrt(1 + tau^2));
    else
        t = -1/(-tau + sqrt(1 + tau^2));
    end
    c = 1/sqrt(1 + t^2);
    s = t*c;
else
    c = 1;
    s = 0;
end