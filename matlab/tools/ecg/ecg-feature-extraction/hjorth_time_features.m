function feature = hjorth_time_features(data, fs)
    % feature = hjorth_time_features(data , fs)    
    % Extract features related to the complexity analysis of ECG signal
    %
    % Inputs:
    %   data: ECG signal (1D array).
    %
    % Output:
    %   feature: Structure containing three values (activity, mobility and complexity)
    %
    % Author:
    %   Seyedeh Somayyeh Mousavi 
    %   Emory University, Georgia, USA
    %   Email: bmemousavi@gmail.com
    %   Date: SEP 26, 2024
    %   Sajjad Karimi
    %   Date: Mar 11, 2025


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

    % Return mobility and complexity
    feature.activity = x_var;
    feature.mobility = mob;
    feature.complexity = com;
end