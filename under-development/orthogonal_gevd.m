function [W, lambda] = orthogonal_gevd(A, B)

if size(A, 1) ~= size(A, 2) || size(B, 1) ~= size(B, 2) || size(A, 1) ~= size(B, 1)
    error('A and B must be square matrices of the same size')
end

n = size(A, 1);

W = zeros(n);
lambda = zeros(n, 1);

A_proj = A;
B_proj = B;
for k = 1 : n
    [v, lambda(k)] = eigs(A_proj, B_proj, 1, 'largestabs');
    lambda(k) = lambda(k) * norm(v);
    v = v/norm(v);
    W(k, :) = v';
    if k < n
        P = eye(n) - v * v';
        P = (P + P') / 2; % for numeric stability only
        A_proj = P * A_proj * P;
        B_proj = P * B_proj * P;
    end
end