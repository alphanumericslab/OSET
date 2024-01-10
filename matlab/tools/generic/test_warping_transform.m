% A test script for warping_transform function
clear; close all; clc

% Example 1:
% Generate warping matrix using quadratic interpolation
in_len = 50;
out_len = 57;
interp_order = 2;
M = warping_transform_col(in_len, out_len, interp_order);

% Example 2:
% Define input and output knot points
% in_knots = [1, 50, 100, 500, 1000];
% out_knots = [1, 70, 140, 700, 1000];
in_knots = [1, 500, 1000];
out_knots = [1, 1100, 2000];

% Generate warping matrix using linear interpolation
interp_order = 2;
M = warping_transform_col(in_knots, out_knots, interp_order);
% The intermediate points are guaranteed to be mapped to one another

% Create an input vector
% input_vector = cumsum(rand(1000, 1));
input_vector = cumsum(randn(1000, 1));

% Warp the input vector using the generated matrix
warped_vector = M * input_vector;

figure
in_time = (0:length(input_vector)-1);
out_time = length(input_vector)*(0:length(warped_vector)-1)/length(warped_vector);

figure
plot(in_time, input_vector);
hold on
plot(out_time, warped_vector);
plot(in_time(in_knots), input_vector(in_knots), 'kx', 'markersize', 18);
grid

% figire
% plot()

