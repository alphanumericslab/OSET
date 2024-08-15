function channel_feature = functioncaller_1(method_name, data, k, file_name, Params)

% functioncaller_1 - Calls a specified function to detect peaks in data.

% Inputs:
%   method_name: Name of the function that you want to call
%   input_data: Data on which method is performed
%   k: Index of channel
%   file_name: Name of the file to load data from
%   data_input_path: Folder where the input files are located
%   Params: Structure containing parameters for feature extraction

% Output:
%   channel_feature: Struct containing peak detection results and success status

% Date: June 28, 2024
% Location: Emory University, Georgia, USA
% By: Seyedeh Somayyeh Mousavi
% Email: bmemousavi@gmail.com 

% Initialize the channel feature struct
channel_feature = struct();
channel_feature.success = 0;
channel_feature.warning = [];

try
    data_channel = data(k, :);
    % Call the function to detect peaks
    [peaks, peak_indexes_consensus] = feval(method_name, data_channel, Params.Fs);
    
    % Store results and indicate success
    channel_feature.success = 1 ;
    channel_feature.peaks = peaks;
    channel_feature.peak_indexes_consensus = peak_indexes_consensus;

catch ME
    % Handle errors during peak detection
    channel_feature.warning = warning(['Failed to perform Func: ' method_name ' for file:' file_name '. Error:' ME.message]);
    
    % Store empty results and indicate failure 
    channel_feature.peaks = [];
    channel_feature.peak_indexes_consensus = [];
end

end
