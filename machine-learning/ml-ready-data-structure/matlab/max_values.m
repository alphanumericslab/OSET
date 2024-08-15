function [maxvalues] = max_values(data)

    % max_values: Function to find the maximum value in each column of the input data matrix

    % Input:
    %   data - A matrix of any size where we want to find the maximum value in each column

    % Output:
    %   maxvalues - A row vector containing the maximum value in each column of the input data

    % Find the maximum values along the first dimension (rows)

    % Date: June 23, 2024
    % Location: Emory University, Georgia, USA
    % By: Seyedeh Somayyeh Mousavi
    % Email: bmemousavi@gmail.com 
    
    maxvalues = max(data, [], 1);
end