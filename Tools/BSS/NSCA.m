function [y, W, A, B, C] = NSCA(x,I,J) % nonstationary component analysis
%
% I, J: desired time windows

B = cov(x(:,I)');
C = cov(x(:,J)');

C = (C+C')/2;
B = (B+B')/2;

[V,D] = eig(B,C,'chol');

d = diag(D);
[YY, II] = sort(d);
II = II(end:-1:1);

W = V(:,II)';
A = pinv(W);

y = real(W*x);

