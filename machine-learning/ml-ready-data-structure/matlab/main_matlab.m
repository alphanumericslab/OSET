% General template of MATLAB script for extracting features using various functions from signals 
% 
% Date: June 19, 2024
% Location: Emory University, Georgia, USA
% By: Seyedeh Somayyeh Mousavi
% Email: bmemousavi@gmail.com

% Clear command window, workspace, and close all figures
clc         % Clear command window
clear       % Clear all variables from workspace
close all   % Close all open figures

%% Part 1: Setup and Initialization (Modify)

% Part 1.1: Path containing the data (Modify)
Data_Input_Path = '../sample_data/';

% Part 1.2: Path where the results will be saved (Modify)
Data_Output_Path = '../results/';

% Part 1.3: Name of the output file (Modify)
Desired_Output_File = 'Example_ECG_ML_Ready_Matlab.json';

% Part 1.4: List of files in the input folder and subfolders (Modify)
Input_Data_list = dir(fullfile(Data_Input_Path, '**/*.mat'));

% Part 1.5: List of MD5_Checksum files in the output folder (Don't change)
Input_Checksum_list = dir(fullfile(Data_Output_Path, 'MD5_Checksum_records_matlab', '**/*.txt'));

%% Part 2: Parameter Configuration (Modify)

% Part 2.1: Sampling frequency (Hz)
Params.Fs = 500.0;

% Part 2.2: Number of channels in the dataset
Params.Number_Channels = 12;

% Part 2.3: Powerline frequency (Hz)
Params.Powerline_frequency = 60;

%% Part 2: Parameter Configuration (Don't change)

% Part 2.4: Define the function details
filePath = fullfile(Data_Output_Path, 'MD5_Checksum_records_matlab', 'Params_function_dict.csv');

% Reading the CSV file into a table
data = readtable(filePath); 

% Initialize Params structure
Params.function_dict = struct();

% Loop through each row in the table and populate Params structure
for i = 1:height(data)
    % Create the function name based on the function_key column
    func_key = data.function_key{i} ;
    Params.function_dict.(func_key) = struct();
    
    %  Assign values directly since they are either strings or empty
    Params.function_dict.(func_key) = struct( ...
        'name', data.name(i), ...
        'rel_path', data.rel_path(i), ...
        'codebase_md5chsum', data.codebase_md5chsum(i), ...
        'codebase_git_repo', data.codebase_git_repo(i), ...
        'codebase_git_commit_id', data.codebase_git_commit_id(i) ...
    );
end

%% Part 3: Validations (Don't Change)

% Check if there are no files in the input data list
if isempty(Input_Data_list)
   error('No signal files found in the Data_Input_Path. Please check the path and try again.');
end

% Check if Fs is a numerical value
if ~isnumeric(Params.Fs)
   error('Fs must be a numerical value.');
end

% Check if Number_Channels is a numerical value and an integer
if ~isnumeric(Params.Number_Channels) || mod(Params.Number_Channels, 1) ~= 0
   error('Number of channels must be an integer.');
end

% Check if the Output_Path exists
if ~isfolder(Data_Output_Path)
   error('Output_Path "%s" does not exist. Please check the path and try again.', Data_Output_Path);
end

%% Part 4: Feature Extraction (Don't Change)

% Call the feature extraction function with the specified parameters
% - Input_Data_list: List of files to process
% - Input_Checksum_list: List of MD5 checksums for the data
% - Params: Struct containing various parameters including sampling
%   frequency, number of channels, and all input variables of the feature extraction functions
% - Data_Output_Path: Folder to save the results
% - Desired_Output_File: Name of the output file

extract_feature(Input_Data_list, Input_Checksum_list, Params, Data_Output_Path, Desired_Output_File);
