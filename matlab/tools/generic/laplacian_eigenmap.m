function [V, d, delta, similarity, epsilon] = laplacian_eigenmap(C, kappa)
% laplacian_eigenmap - Calculates Laplacian eigen-maps of positive semi-definite matrix pairs.
%
% Syntax: [V, d, delta, similarity, epsilon] = laplacian_eigenmap(C, kappa)
%
% Inputs:
%   C: K semi-positive definite matrices of size NxN provided as a N*N*K tensor.
%   kappa: The median factor in distance measure (see implementation).
%
% Outputs:
%   V: Eigenvectors of the Laplacian eigen-maps.
%   d: Corresponding eigenvalues of the Laplacian eigen-maps.
%   delta: Eigen-map distance matrix.
%   similarity: Eigen-map similarity matrix.
%   epsilon: kappa times the median of delta.
%
%   Revision History:generalized
%       2020: First release
%       2023: Renamed from deprecated version LaplacianEigenmap()
%
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

K = size(C, 3);   % Number of matrices in the tensor
delta = zeros(K);   % Initialize the eigen-map distance matrix

% Compute pairwise eigen-map distance
for i = 1 : K-1
    for j = i + 1 : K
        lambdas = eig(C(:, :, i), C(:, :, j));   % Compute eigenvalues
        delta(i, j) = sum(log(abs(lambdas)).^2);   % Compute eigen-map distance
        delta(j, i) = delta(i, j);   % Symmetric matrix
    end
end

% Compute epsilon and similarity matrix
mask = logical(triu(ones(K), 1));   % Upper triangular mask
dd = delta(mask);   % Extract upper triangular elements
epsilon = kappa * median(dd(:));   % Compute epsilon
similarity = exp(-delta / epsilon);   % Compute similarity matrix

% Compute Laplacian eigen-maps
S = diag(1 ./ sqrt(sum(similarity, 2)));   % Compute diagonal scaling matrix
[V, D] = eig(S * similarity * S);   % Compute eigenvectors and eigenvalues
[d, I] = sort(diag(D), 'descend');   % Sort eigenvalues in descending order
V = V(:, I);   % Reorder eigenvectors

