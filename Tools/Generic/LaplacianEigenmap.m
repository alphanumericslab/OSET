function [V, d, delta, similarity, epsilon] = LaplacianEigenmap(C, kappa)
% [V, d, delta, similarity, epsilon] = LaplacianEigenmap(C, kappa)
% Calculates Laplacian eigen-maps of positive semi definite matrix pairs
%
% Inputs:
%   C: K semi-positive definite matrices of size NxN provided as a N*N*K tensor
%   kappa: The median factor in distance measure (see implementation)
% 
% Outputs:
%   V: eigenvectors of the Laplacian eigen-maps
%   d: corresponding eigenvalues of the Laplacian eigen-maps
%   delta: eigen-map distance matrix
%   similarity: eigen-map similarity matrix
%   epsilon: kappa*median(delta) 
% 
% Reza Sameni, 2020
% 
% The The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET.git
% 
% History: Updated Jan, 2023


K = size(C, 3);
delta = zeros(K);
for i = 1 : K-1
    for j = i + 1 : K
        lambdas = eig(C(:, :, i), C(:, :, j));
        delta(i, j) = sum( log(abs(lambdas)) .^ 2);
        delta(j, i) = delta(i, j);
    end
end
% delta = delta + delta';
mask = logical(triu(ones(K), 1));
dd = delta(mask);
epsilon = kappa * median(dd(:));
similarity = exp(-delta / epsilon);
S = diag(1 ./ sqrt(sum(similarity, 2)));
[V, D] = eig(S*similarity*S);
[d, I] = sort(diag(D), 'descend');
V = V(:, I);