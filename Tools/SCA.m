function [y,W,A] = SCA(x,f1,f2) % spectral component analysis
%
% f1,f2: Normalized frequencies

L1 = size(x,1);
L2 = size(x,2);

k = max(round(L2*f1), 1) : min(round(L2*f2), L2);

k = [k L2-k+1];
Ax = zeros(L1,L1);
for ii = 1:L1
    X = fft(x(ii,:),L2);
    for jj = 1:L1
        Y = fft(x(jj,:),L2);
        c = X.*conj(Y);
        Ax(ii,jj) = sum(c(k));
    end
end
B = real(Ax);
C = x*x';

C = (C+C')/2;
B = (B+B')/2;

[V,D] = eig(B,C,'chol');

d = diag(D);
[YY,II] = sort(d);
II = II(end:-1:1);

W = V(:,II)';
A = pinv(W);

y = real(W*x);

