function [y,W,A] = SCA2(x, f0, bw, order) % spectral component analysis based on BP filtering
%
% f1,f2: Normalized frequencies

% xx = BPFilter(x, f1, f2);
xx = BPFilter5(x, f0, bw, order);

% load bandpassfilter01 h
% xx = filter(h, 1, x, [], 2);


B = xx*xx';
C = x*x';

B = (B+B')/2;
C = (C+C')/2;

[V,D] = eig(B,C,'chol');

d = diag(D);
[~,II] = sort(d);
II = II(end:-1:1);

W = V(:,II)';
A = pinv(W);

y = real(W*x);

