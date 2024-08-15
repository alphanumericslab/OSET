# MATLAB-Based Feature Extraction for AI and ML Learning Applications

## Contents
This folder contains MATLAB scripts for applying various feature extraction functions to a dataset. Additionally, using a simple example, we explain how to modify the scripts based on your specific needs.

## Run the Codes
1. Clone the [OSET repository](https://github.com/alphanumericslab/OSET).
2. Add the main folder and its subfolders to your MATLAB path.
3. Run [main_matlab](main_matlab.m). You will see the JSON file in the output folder.
4. The following table lists the main MATLAB scripts of the feature extraction pipeline. Modify these scripts to run on your dataset based on the feature extraction functions you have considered.

| File Name                   | Description                                                                                                    | 
|-----------------------------|----------------------------------------------------------------------------------------------------------------|
| [main_matlab](main_matlab.m)  | The script initializes parameters, searches for data files, validates inputs, and calls a feature extraction function to process the data, then saves the results.| 
|[extract_feature](extract_feature.m) | The script loads the data and then runs the feature extraction algorithms on it.   | 
|functioncaller |The script includes a function that calls a specific function.|

## Modify the Codes
1. In the first step, before running the MATLAB code, you should generate the MD5 checksums of the data and functions that you have considered. Open the [Generate_MD5_Checksum](Generate_MD5_Checksum.py).
2. You should modify `Part 1: Setup and Initialization` to include the path containing the data and the path where the results will be saved. Additionally, you should define the type of data, such as .mat files.
3. Then, you should modify `Part 3: Parameter Configuration` to define the list of functions needed for feature extraction. `load` is a necessary function (Don't change it). Please add the other functions in the format provided below to the `Params['function_dict']`. All function inputs should be strings.
   
```
# Change n, the number that you assign to a function (e.g., 1, 2, 3, ...)
Fun_n: {
    'name': 'Fill this part, method/function name (string)', ...
    'rel_path': 'Fill this part, method/function relative path to codebase root (string)', ...
    'codebase_git_repo': 'Fill this part, Git repository URL/link (string)', ...
    'codebase_git_commit_id': 'Fill this part, Unique git commit ID (HEX string)' ...
}
```
**Please ensure that the function names are accurate. It is crucial that they are called exactly as specified.**

`Fun_1` and `Fun_2` are functions related to example codes. You can deactivate or change them.

4. Run the script to generate the MD5 checksums for the data and functions.      
5. Open the [main_matlab](main_matlab.m).
6. You should modify `Part 1: Setup and Initialization`. According to the aim of your project, you should define the following 4 parts:
   - Data_Input_Path: Path containing the data
   - Data_Output_Path: Path where the results will be saved
   - Desired_Output_File: Name of the output file
   - The given extension (for example, .mat, .csv, or .dat) to list the files in the folder and subfolders of Data_Input_Folder
7. Then, you should modify `Part 2: Parameter Configuration` to correct the sampling frequency, number of channels, and powerline frequency. 
8. Open the [extract_feature](extract_feature.m).
9. You need to write scripts to call the functions and then save the outputs. Each function has four parts:
- Part 1: Function Setup and Initialization
```matlab
% Change n, the number related to the function that you have defined in main_matlab file . 
method = list_functions.Fun_n;
% Don't change
Num_Fun = Num_Fun + 1;
feature.methods.(['method_' num2str(N_Fun)]).name = method.name;
feature.methods.(['method_' num2str(N_Fun)]).rel_path = method.rel_path;
feature.methods.(['method_' num2str(N_Fun)]).codebase_md5chsum = method.codebase_md5chsum;
feature.methods.(['method_' num2str(N_Fun)]).codebase_git_repo = method.codebase_git_repo;
feature.methods.(['method_' num2str(N_Fun)]).codebase_git_commit_id = method.codebase_git_commit_id;
```
- Part 2: Call the function, choose the appropriate option based on whether you want to use whole data or process channels one by one (Whole_data, Channel_data)

```matlab
% Part 2: Whole_data
% Change n, the number related to the function that you have defined in main_matlab file .  
result = functioncaller_n(data, file_name, N_Fun, Params);

% Store the input info of the function
% Don't change
feature.methods.(['method_' num2str(N_Fun)]).inputs.name = 'Whole_data';
% Change 
feature.methods.(['method_' num2str(N_Fun)]).inputs.type = [Fill this part, input type (string): 'bool', 'int', 'float' or 'string'];
% Don't change. We considered 0 for the list of all channels from the input
feature.methods.(['method_' num2str(N_Fun)]).inputs.channels = 0;
```
```matlab
% Part 2: Channel_data
% Don't change 
channel_name = cell(1, Params.Number_Channels);
warning = cell(1, Params.Number_Channels);
success = zeros(1, Params.Number_Channels);

% Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
% ** output = cell(1, Params.Number_Channels) **
output1 = cell(1, Params.Number_Channels);
output2 = cell(1, Params.Number_Channels);
        
% Process each channel individually
for k = 1:Params.Number_Channels
   % Call the function
   % Change n, the number related to the function that you have defined in main_matlab file .
   result_n = functionCaller_n(data, k, file_name, N_Fun, Params);
   % Don't chanage  
   channel_name{k} = k;
   % Change the n, number related to the function that you have defined in main_matlab file .
   success(k) = result_n.success;
   warning{k} = result_n.warning;
   % Change n 
   % Add the line of code below, based on the number of the outputs, and replace the placeholders with the internal names of the outputs
   % ** output{k} = result_n.output **
   output1{k} = result_n.output1;
   output2{k} = result_n.output2;
end

% Store the input info of the function
% Don't change
feature.methods.(['method_' num2str(N_Fun)]).inputs.name = 'Channel_data';
% Change 
feature.methods.(['method_' num2str(N_Fun)]).inputs.type = [Fill this part, input type (string): 'bool', 'int', 'float' or 'string'];
% Don't change
feature.methods.(['method_' num2str(N_Fun)]).inputs.channels = channel_name;
```
- Part 3: Save the parameters.
```
% Part 3: Store the input parameters of the function
% You should modify the line codes and replace the placeholders with the internal name of the parameter
feature.methods.(['method_' num2str(N_Fun)]).params.name = []; % Use internal name of the parameter
feature.methods.(['method_' num2str(N_Fun)]).params.type = []; % 'bool', 'int', 'float' or 'string'
feature.methods.(['method_' num2str(N_Fun)]).params.val = []; % bool/int/float/string, single-valued or array
```
If there are multiple parameters, replace `params` with `params.param_N`, where N represents the parameter number (e.g., 1, 2, ...). Repeat this structure for each additional parameter.

- Part 4: Save the results, choose the appropriate option based on Part 2 (Whole_data, Channel_data)
If there are multiple outputs, replace `outputs` with `outputs.output_N`, where N represents the parameter number (e.g., 1, 2, ...). Repeat the structure for each output.

```matlab
% Part 4: Whole_data
% Save the result of the function
% You should modify the line codes.
feature.methods.(['method_' num2str(N_Fun)]).outputs.name = []; % Use internal variable name of the output
feature.methods.(['method_' num2str(N_Fun)]).outputs.type = []; % output type (string): 'bool', 'int', 'float' or 'string'
% Change n and replace the placeholder with the internal name of the output
feature.methods.(['method_' num2str(N_Fun)]).outputs.val = result_n.output;
% Change n 
feature.methods.(['method_' num2str(Num_Fun)]).warning = result_n.warning;
feature.methods.(['method_' num2str(Num_Fun)]).success = result_n.success;
```
```
% Part 4: Channel_data
% Save the result of the function
% You should modify the line codes.
feature.methods.(['method_' num2str(N_Fun)]).outputs.name = []; % Use internal variable name of the output
feature.methods.(['method_' num2str(N_Fun)]).outputs.type = []; % output type (string): 'bool', 'int', 'float' or 'string'
% Replace the placeholder with the internal name of the output
feature.methods.(['method_' num2str(N_Fun)]).outputs.val = output;
% Don't change 
feature.methods.(['method_' num2str(Num_Fun)]).warning = warning;
feature.methods.(['method_' num2str(Num_Fun)]).success = success;
```
10. For each feature extraction function, copy the four parts and modify them according to the function's name and its outputs.
11. Each function has a `functioncaller` file that includes a function to call a specific function. Choose the script based on whether you want to use the Whole_data or Channel_data, and create the corresponding functioncaller_n file where n represents the number related to the function that you have defined in main_matlab file.  
    
**Given the importance of preventing code crashes, it's crucial to have a good familiarity with the input parameters, outputs, and types of them. Therefore, ensure to include both parts of the try and catch blocks in the code.**
    
```matlab
% Whole_data
% Change n, the number related to the function that you have defined in main_matlab file . 
function feature = functioncaller_n(method_name, data, file_name, Params)

% Initialize the output structure (Don't change)
feature = struct();
feature.warning = [];
feature.success = 0;

try
    % Modify Params.var and output based on the function you want to call. Use the internal names of the parameters and outputs.
    [output_1, output_2, ...] = feval(method_name, data, Params.var_1, Params.var_2, ...);
    
    % Save results and indicate success (Don't change)
    feature.success = 1;
    % Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
    % Modify the output based on the function you want to call. Use the internal names of the outputs.
    % ** feature.output_N = output_N; **
    feature.output_1 = output_1;
    feature.output_2 = output_2;
    
catch ME
    % Handle any errors that occur during function execution (Don'T change)
    feature.warning = warning(['Failed to perform Func: ' method ' for file: ' file_name '. Error: ' ME.message]);
    
    % Store empty results and indicate failure
    % Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
    % Modify the output_N based on the function you want to call. Use the internal names of the outputs.
    % ** feature.output_N = []; **
    feature.output_1 = [];
    feature.output_2 = [];
end
end
```

```matlab
% Channel_data
% k is the index of channel
% Change n, the number related to the function that you have defined in main_matlab file .  
function channel_feature = functioncaller_n(method_name, data, k, file_name, Params)

% Initialize the channel feature struct (Don't change)
channel_feature = struct();
channel_feature.success = 0;
channel_feature.warning = [];

try
    % k is the index of channel  
    data_channel = data(k, :);
    % Modify Params.var and output based on the function you want to call. Use the internal names of the parameters and outputs.
    [output_1, output_2, ...] = feval(Params.list_functions{N_Fun}, data_channel, Params.var_1, Params.var_2, ...);
    
    % Save results and indicate success (Don't change)
    channel_feature.success = 1;
    % Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
    % Modify the output based on the function you want to call. Use the internal names of the outputs.
    % ** channel_feature.output_N = output_N; **
    channel_feature.output_1 = output_1;
    channel_feature.output_2 = output_2;
    
catch ME
    % Handle errors during peak detection (Don't change)
    channel_feature.warning = warning(['Failed to perform Func: ' method ' for file:' file_name '. Error:' ME.message]);
 
    % Save empty results and indicate failure (Don't change)
    % Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
    % Modify the output based on the function you want to call. Use the internal names of the outputs.
    % ** channel_feature.output_N = []; ** 
    channel_feature.output_1 = [];
    channel_feature.output_2 = [];
end
end
```  
12. `load` function is essential. Open [functioncaller_0](https://github.com/alphanumericslab/ecg-blood-pressure/blob/main/ml-ready-data-structure/matlab/functioncaller_0.m) and modify the code to correctly read the data. Each dataset has a specific format, so ensure the loading process is adapted to handle your dataset's format.
13. After modifying the code and creating the relevant `functioncaller_n` files, you can easily run [main_matlab](https://github.com/alphanumericslab/ecg-blood-pressure/blob/main/ml-ready-data-structure/matlab/main_matlab.m) in the cluster, AWS, or any other environments. You will see the JSON file in the output folder. 

## Contributors
- [Somayyeh Mousavi](seyedeh.somayyeh.mousavi@emory.edu)
- [Reza Sameni](rsameni@dbmi.emory.edu)