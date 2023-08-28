function R = Rotate3D(theta_x, theta_y, theta_z)
% Rotate3D has been deprecated. Use rotation_matrix_3d instead.
warning('Rotate3D has been deprecated. Use rotation_matrix_3d instead.');
R = rotation_matrix_3d(theta_x, theta_y, theta_z);