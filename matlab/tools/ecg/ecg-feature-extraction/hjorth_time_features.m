function [feature_vec, feature_info] = hjorth_time_features(data, fs)

% [feature_vec, feature_info] = hjorth_time_features(data)
% Extract features related to the complexity analysis of ECG signal
%
% Inputs:
%   data: ECG signal as a vector (in microvolts) (1D array).
%
% Output:
%   feature_vec:  A vector contains three features (activity, mobility and complexity)
%   feature_info: A structure contains feature descriptions, names and
%   units.
%
% Author:
%   Seyedeh Somayyeh Mousavi
%   Sajjad Karimi
%   Reza Sameni
%   Emory University, Georgia, USA
%   Email: bmemousavi@gmail.com
%   First: Date: SEP 24, 2024
%   Second: Date: AUG 3, 2025

% Calculate derivatives
dx = fs*diff(data);     % First derivative
ddx = fs*diff(dx);      % Second derivative

% Calculate variance
x_var = var(data);      % Variance of original signal
dx_var = var(dx);       % Variance of first derivative
ddx_var = var(ddx);     % Variance of second derivative

% Mobility and complexity calculations
mob = sqrt(dx_var / x_var);             % Mobility
com = sqrt(ddx_var / dx_var) / mob;     % Complexity

% Return activity, mobility and complexity
% feature.activity = x_var;
% feature.mobility = mob;
% feature.complexity = com;

feature_vec = [x_var, mob, com];
% Define feature info
feature_info.names = {'activity','mobility', 'complexity'};
feature_info.units = {'mv^2','Hz', 'scaler'};
feature_info.description = {"Hjorth parameters: Activity","Hjorth parameters: Mobility","Hjorth parameters: Complexity"};

end