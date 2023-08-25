function R = rotation_matrix_3d(theta_x, theta_y, theta_z)
% rotation_matrix_3d - Generate a 3D rotation matrix.
%
% Syntax:
%   R = rotation_matrix_3d(theta_x, theta_y, theta_z)
%
% Inputs:
%   theta_x: Rotation angle around the x-axis or in y-z plane (in radians)
%   theta_y: Rotation angle around the y-axis or in x-z plane (in radians)
%   theta_z: Rotation angle around the z-axis or in x-y plane (in radians)
%
% Output:
%   R: The rotation matrix as the product of three rotation matrices around
%      x, y, and z (in order): R = Rx * Ry * Rz
%
% Description:
%   This function generates a 3D rotation matrix by performing rotations
%   around the x, y, and z axes. The resulting matrix represents the composite
%   rotation.
%
% Revision History:
%   2006: First release.
%   2023: Documented and renamed from deprecated version Rotate3D
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Generate rotation matrices around the x, y, and z axes.
Rx = [1, 0, 0; 0, cos(theta_x), sin(theta_x); 0, -sin(theta_x), cos(theta_x)];
Ry = [cos(theta_y), 0, -sin(theta_y); 0, 1, 0; sin(theta_y), 0, cos(theta_y)];
Rz = [cos(theta_z), sin(theta_z), 0; -sin(theta_z), cos(theta_z), 0; 0, 0, 1];

% Calculate the composite rotation matrix.
R = Rx * Ry * Rz;
