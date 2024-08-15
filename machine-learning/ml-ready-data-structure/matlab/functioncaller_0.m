function [feature, data] = functioncaller_0(method_name, file_name, data_input_path, Params)

% functioncaller_0 - Calls load function

% Inputs:
%   method_name: Name of the function that you want to call
%   file_name: Name of the file to load data from
%   data_input_path: Folder where the input files are located
%   Params: Structure containing parameters for feature extraction


% Outputs:
%   feature: Struct indicating success of the function call (success = true/false)
%   data: Loaded data from the file after processing (if successful)

% Date: June 28, 2024
% Location: Emory University, Georgia, USA
% By: Seyedeh Somayyeh Mousavi
% Email: bmemousavi@gmail.com 

% Don't change
feature.warning = [];
data =[];

try
    % Load the data using the specified function (Modify this part)
    data = load(fullfile(data_input_path, file_name));
    
    % Modify this part: Extract the 'val' field from loaded data
    data = data.val; 
    
    % Don't change
    feature.success = 1;

catch ME
    % Don't change
    % Handle errors during data loading 
    feature.warning= warning(['Failed to perform Func: ' method_name ' for file: ' file_name '. Error: ' ME.message]);
    feature.success = 0; 
end

% Don't change
% Process data if loading was successful
if feature.success
    % Check the shape of the data and reshape if necessary
    if size(data, 1) ~= Params.Number_Channels
        data = data'; % Transpose data if number of rows doesn't match Number_Channels
    end
end

end
