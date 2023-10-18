function [W_a, lambda_a, W_b, lambda_b] = orthogonal_gevd(A, B)

if size(A, 1) ~= size(A, 2) || size(B, 1) ~= size(B, 2) || size(A, 1) ~= size(B, 1)
    error('A and B must be square matrices of the same size')
end

n = size(A, 1);

W_a = zeros(n);
lambda_a = zeros(n, 1);

W_b = zeros(n);
lambda_b = zeros(n, 1);


A_proj = A;
B_proj = B;
for k = 1 : n
    [v, lambda_a(k)] = eigs(A_proj, B_proj, 1, 'largestabs');
    lambda_a(k) = lambda_a(k) * norm(v);
    v = v/norm(v);
    W_a(k, :) = v';
    if k < n
        P = eye(n) - v * v';
        P = (P + P') / 2; % for numeric stability only
        A_proj = P * A_proj * P;
        B_proj = P * B_proj * P;
    end
end

A_proj = A;
B_proj = B;
for k = 1 : n
    [v, lambda_b(k)] = eigs(B_proj, A_proj, 1, 'largestabs');
    lambda_b(k) = lambda_b(k) * norm(v);
    v = v/norm(v);
    W_b(k, :) = v';
    if k < n
        P = eye(n) - v * v';
        P = (P + P') / 2; % for numeric stability only
        A_proj = P * A_proj * P;
        B_proj = P * B_proj * P;
    end
end