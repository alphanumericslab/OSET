function [P_tilde, P] = mica_projectors(A)
% mica_projectors - Multidimensional Independent Component Analysis (MICA) projection matrices.
% 
% Description:
%   Given an nxn matrix A, which has been estimated by, for example,
%   independent component analysis of mixtures following the model x = A*s,
%   this function provides the projectors that can be used to decompose x
%   into its linear decomposition as follows: x = \sigma_p x_p. Using the
%   projectors we can find: x_p = P_tilde{p} * x. See the reference below 
%   for details and notations.
% 
% Usage: [P_tilde, P] = mica_projectors(A);
% 
% Inputs:
%   A: An nxn matrix representing the estimated mixing matrix in the ICA model.
%
% Outputs:
%   P_tilde: A cell array of n projection matrices (each n x n) for each independent component.
%   P: A cell array of n projection matrices (each n x n) used to create P_tilde.
% 
% Reference:
%   Cardoso, J-F. "Multidimensional independent component analysis." In
%   Proceedings of the 1998 IEEE International Conference on Acoustics,
%   Speech and Signal Processing, ICASSP'98 (Cat. No. 98CH36181), vol. 4,
%   pp. 1941-1944. IEEE, 1998. doi: 10.1109/ICASSP.1998.681443
% 
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

n = size(A, 2); % Get the number of columns in A (assumes square matrix).

P = cell(1, n); % Initialize a cell array for P.
S = zeros(n);   % Initialize S as an n x n matrix of zeros.

for k = 1 : n
   % Calculate the projection matrix P{k} for each independent component.
   P{k} = A(:, k) * A(:, k)' / (A(:, k)' * A(:, k));
   S = S + P{k}; % Accumulate P{k} into S.
end

P_tilde = cell(1, n); % Initialize a cell array for P_tilde.
S_inv = pinv(S);      % Calculate the pseudo-inverse of S.

for k = 1 : n
    % Calculate P_tilde by multiplying each P{k} with S_inv.
    P_tilde{k} = P{k} * S_inv;
end

end
