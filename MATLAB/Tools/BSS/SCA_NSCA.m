function [y, W, A] = SCA_NSCA(x,f1,f2, I, order) % sort components by spectrally desired components vs nonstationary components
%
% f1,f2: Normalized frequencies
% I: desired/underired indexes
% order: 'forward' vs 'backward' sorting of spectral and nonstationary components

% The covariance over the frequency band of interest
xx = BPFilter(x, f1, f2);
B = xx*xx';
B = (B+B')/2;

% The covariance over the indexes of interest
C = cov(x(:, I)');
C = (C+C')/2;

if isequal(order, 'forward')
    [V,D] = eig(B, C, 'chol');
elseif isequal(order, 'backward')
    [V,D] = eig(C, B, 'chol');
end

d = diag(D);
[~,II] = sort(d);
II = II(end:-1:1);

W = V(:,II)';
A = pinv(W);

y = real(W*x);

