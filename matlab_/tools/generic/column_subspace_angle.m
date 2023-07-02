function [ang, theta_min, theta_max, theta_mean, theta_med] = column_subspace_angle(A, B)
% column_subspace_angle - Compute angle between the column subspaces of matrices A and B.
%
% Syntax: [ang, theta_min, theta_max, theta_mean, theta_med] = column_subspace_angle(A, B)
%
% Inputs:
%   A: Matrix A with dimensions M x N1.
%   B: Matrix B with dimensions M x N2.
%
% Outputs:
%   ang: The column subspace angles (N1 x N2).
%   theta_min: Minimum angle between column subspaces of A and B.
%   theta_max: Maximum angle between column subspaces of A and B.
%   theta_mean: Mean angle between column subspaces of A and B.
%   theta_med: Median angle between column subspaces of A and B.
%
%   Revision History:
%       2015: First release
%       2023: Renamed from deprecated version ColumnSubspaceAngle()
%
%   Reza Sameni, 2015-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if(size(A, 1) ~= size(B, 1))
   error('A and B should have the same number of rows'); 
end

ang = zeros(size(A, 2), size(B, 2));
for i = 1 : size(A, 2)
    for j = 1 : size(B, 2)
        ang(i, j) = 180*acos(A(:,i)'*B(:,j)/sqrt((A(:,i)'*A(:,i)) * (B(:,j)'*B(:,j)))) / pi;  % Compute the angle between column subspaces
    end
end

theta_min = min(ang(:));    % Compute the minimum angle
theta_max = max(ang(:));    % Compute the maximum angle
theta_mean = mean(ang(:));  % Compute the mean angle
theta_med = median(ang(:)); % Compute the median angle
