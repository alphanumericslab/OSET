function [S, R, e] = subspace_mapping(A, B, ITR, mode)
% SUBSPACE_MAPPING  Estimates scaling (S) and rotation (R) matrices between subspaces A and B.
%
%   [S, R, e] = subspace_mapping(A, B, ITR) finds diagonal scaling S and
%   orthogonal rotation R that minimize the Frobenius norm:
%
%       ||A - B*S*R||_F
%
%   Inputs:
%       A, B : Matrices of equal size (m×n) representing subspaces.
%       ITR  : Number of iterations for better convergence (default = 1).
%       mode : 'rotate' (det(R) = +1), 'reflect' (det(R) = -1) or 'any' (det(R) = +/-1) 
%
%   Outputs:
%       S : Diagonal scaling matrix (n×n)
%       R : Rotation matrix (n×n)
%       e : Final Frobenius norm error
%
%   Reference: Inspired by the subspace rotation method in:
%       Golub, G.H. & Van Loan, C.F. (1996). Matrix Computations, Chapter 12.4.
%
%   Reza Sameni, 2025
%   The Open-Source Electrophysiological Toolbox (OSET)
%   https://github.com/alphanumericslab/OSET

if size(A,1) ~= size(B,1) || size(A,2) ~= size(B,2)
    error('A and B must have the same size');
end
if nargin < 3 || isempty(ITR)
    ITR = 1;
end

if nargin < 4 || isempty(mode)
    mode = 'any';
end


n = size(A,2);
S = eye(n);

for k = 1:ITR
    [U, sigma, V] = svd(S * B' * A);
    R = U * V';
    switch mode
        case 'rotate'
            if det(R) < 0 % flip the sign of the eigenvector corresponding to the smallest singular value
                [~, I] = min(diag(sigma));
                V(:, I) = -V(:, I);
                R = U * V';
            end
        case 'reflect'
            if det(R) > 0 % flip the sign of the eigenvector corresponding to the smallest singular value
                [~, I] = min(diag(sigma));
                V(:, I) = -V(:, I);
                R = U * V';
            end
        case 'any'

        otherwise
            error('Undefined mode')
    end
    Num = B' * A * R';
    Den = B' * B;
    S = diag(diag(Num) ./ diag(Den));
end

e = norm(A - B * S * R, 'fro');
