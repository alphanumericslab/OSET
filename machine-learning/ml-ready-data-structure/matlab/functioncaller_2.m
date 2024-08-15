function feature = functioncaller_2(method_name, input_data, file_name, Params)


% functioncaller_2 - Calls a specified function to calculate maximum values from data.

% Inputs:
%   method_name: Name of the function that you want to call
%   input_data: Data on which method is performed
%   file_name: Name of the file to load data from
%   data_input_path: Folder where the input files are located
%   Params: Structure containing parameters for feature extraction

% Output:
%   feature: Struct containing maximum values and success status

% Date: June 28, 2024
% Location: Emory University, Georgia, USA
% By: Seyedeh Somayyeh Mousavi
% Email: bmemousavi@gmail.com 

% Initialize the output structure
feature = struct();
feature.warning = [];
feature.success = 0;

try
    % call the function to calculate maximum values from data.
    maxvalues = feval(method_name, input_data);
    
    % Store results and indicate success
    feature.success = 1;
    feature.maxvalues = maxvalues;
    
catch ME
    % Handle any errors that occur during function execution
    feature.warning = warning(['Failed to perform Func: ' method_name ' for file: ' file_name '. Error: ' ME.message]);
    
    % Store empty results and indicate failure
    feature.maxvalues = [];
end

end

