% Implementations of the nearest orthogonal matrix (the Orthogonal
% Procrustes problem)
% Refer to: https://en.wikipedia.org/wiki/Orthogonal_matrix
% 
% Reza Sameni, Jan 2019

N = 5
M = randn(N)

% The SVD approach
[U,S,V] = svd(M);
Q = U*V'

% The SVD approach
[U,S,V] = svd(M);
R1 = U*V'

% The square root method
R2 = M*(M'*M)^(-0.5)

% The recursive method
Q = M;
for k = 1 : 10, 
    Q = 2*M*pinv(pinv(Q)* M + M'*Q);
end
R3 = Q
