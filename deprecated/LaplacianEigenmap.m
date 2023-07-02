function [V, d, delta, similarity, epsilon] = LaplacianEigenmap(C, kappa)
% LaplacianEigenmap has been deprecated. Use laplacian_eigenmap instead.
warning('LaplacianEigenmap has been deprecated. Use laplacian_eigenmap instead.');
[V, d, delta, similarity, epsilon] = laplacian_eigenmap(C, kappa);
