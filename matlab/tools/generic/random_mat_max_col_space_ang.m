function [A, B, ang] = random_mat_max_col_space_ang(angle, dim1, dim2, flag, varargin)
% random_mat_max_col_space_ang - Generate random matrices with a maximum angle
% between their column subspaces.
%
% Syntax:
%   [A, B, ang] = random_mat_max_col_space_ang(angle, dim1, dim2, flag)
%   [A, B, ang] = random_mat_max_col_space_ang(angle, dim1, dim2, flag, Kmax)
%
% Inputs:
%   angle: Maximum angle (in degrees) between column subspaces.
%   dim1: First dimension of matrices.
%   dim2: Second dimension of matrices.
%   flag: Flag indicating the generation technique.
%       - flag = 0: Backward compatibility mode (scales the first dim1-3
%           dimensions and rotates the last 3)
%       - flag = 1: General mode.
%   Kmax (optional): Maximum number of iterations before giving up (default: 10000)
%
% Outputs:
%   A: First matrix (dim1 x dim2)
%   B: Second matrix (dim1 x dim2)
%   ang: Maximum angle between column subspaces
%
% Note: The function tries Kmax times and gives up if the requested angle is not found
% 
% Description:
%   This function generates random matrices A and B in n-dimensional space
%   with a specified maximum angle between their column subspaces
%
% Revision History:
%   2008: First release.
%   2023: Documented and renamed from deprecated version RandomMatrices
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Number of iterations before giving up
if nargin > 4 && ~isempty(varargin{1})
    Kmax = varargin{1};
else
    Kmax = 10000; 
end

k = 0;
if flag == 0 % For backward compatibility; scales the first dim1-3 dimensions and rotates the last 3
    ang = [90, 90, 90];
    while max(ang) > angle && k < Kmax
        A = randn(dim1, 3);
        phi_x = pi * angle * rand / 180;
        phi_y = pi * angle * rand / 180;
        phi_z = pi * angle * rand / 180;
        R = [diag(rand(dim1-3, 1)) * eye(dim1-3, dim1); zeros(3, dim1-3) rotation_matrix_3d(phi_x, phi_y, phi_z)];

        B = R * A;

        [Qa, ~] = qr(A, 0);
        [Qb, ~] = qr(B, 0);
        C = Qa' * Qb;
        [~, S, ~] = svd(C);
        S = min(diag(S), 1);
        S = max(S, -1);
        ang = acos(S) * 180 / pi;

        k = k + 1;
    end

elseif flag == 1 % General technique
    ang = [0, 0, 0];
    while min(ang) < angle && k < Kmax
        A = randn(dim1, dim2);
        B = randn(dim1, dim2);

        [Qa, ~] = qr(A, 0);
        [Qb, ~] = qr(B, 0);
        C = Qa' * Qb;
        [~, S, ~] = svd(C);
        S = min(diag(S), 1);
        S = max(S, -1);
        ang = acos(S) * 180 / pi;

        k = k + 1;
    end

end

if k >= Kmax
    A = [];
    B = [];
    ang = [];
    warning(['Giving up. Requested angle not within ', num2str(Kmax), ' iterations.']);
end
