function mean_values = mean_percentile(data, lower_percentile, upper_percentile)
    % Function to calculate the mean of the data within the specified percentile range for each column
    % Input: 
    %   data - a matrix where each column represents a separate dataset
    %   lower_percentile - the lower percentile (e.g., 2.5 for the 2.5th percentile)
    %   upper_percentile - the upper percentile (e.g., 97.5 for the 97.5th percentile)
    % Output:
    %   mean_values - a row vector containing the mean of the data within the specified percentile range for each column
    %
    % Revision History:
    %   2024: First release
    %
    % Reza Sameni, 2024
    % The Open-Source Electrophysiological Toolbox
    % https://github.com/alphanumericslab/OSET

    if isempty(data)
        error('Input data cannot be empty');
    end

    if lower_percentile < 0 || upper_percentile > 100 || lower_percentile >= upper_percentile
        error('Invalid percentile range. Ensure 0 <= lower_percentile < upper_percentile <= 100');
    end

    % Preallocate the output vector
    [~, cols] = size(data);
    mean_values = zeros(1, cols);
    
    % Loop through each column to calculate the mean within the specified percentile range
    for i = 1:cols
        column_data = data(:, i);  % Extract the i-th column
        
        % Calculate the prc for the current column
        prc = prctile(column_data, [lower_percentile, upper_percentile]);

        % Extract the lower and upper percentile values
        lower_value = prc(1);
        upper_value = prc(2);

        % Select data within the specified percentile range
        selected_data = column_data(column_data >= lower_value & column_data <= upper_value);

        % Calculate the mean of the selected data and store it in the output vector
        mean_values(i) = mean(selected_data);
    end
end
