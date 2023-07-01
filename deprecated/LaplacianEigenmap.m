function [V, d, delta, similarity, epsilon] = LaplacianEigenmap(C, kappa)
% Deprecated: LaplacianEigenmap is deprecated. Use laplacian_eigenmap instead.
warning('Deprecated: LaplacianEigenmap is deprecated. Use laplacian_eigenmap instead.');
[V, d, delta, similarity, epsilon] = laplacian_eigenmap(C, kappa);
