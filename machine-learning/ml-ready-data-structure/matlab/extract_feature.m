function [] = extract_feature(Input_Data_list, Input_Checksum_list, Params, Data_Output_Path, Desired_Output_File)
    % extract_feature: MATLAB function for extracting features from signals
    % to create a structure containing all extracted features and relevant
    % info regrading Input_Data and all functions
    
    % Input arguments:
    %   - Input_Data_list: List of files in the input data folder and subfolders
    %   - Input_Checksum_list: List of MD5 checksums for the data
    %   - Params: Structure containing parameters for feature extraction
    %   - Data_Output_Path: Folder where the output data will be saved
    %   - Desired_Output_File: Name of the desired output file
    
    % Date: June 24, 2024
    % Location: Emory University, Georgia, USA
    % By: Seyedeh Somayyeh Mousavi
    % Email: bmemousavi@gmail.com    

    % Initialize the output variable as a struct (Dont'change)
    analysis = struct();

    % Define the list of methods/functions to apply (Dont'change)
    list_functions = Params.function_dict;

    % Loop through each file and process the data (Dont'change)
    for i = 1:length(Input_Data_list)
        
        % Get file name and folder (Dont'change)
        file_name = Input_Data_list(i).name;
        disp(file_name)
        file_path = Input_Data_list(i).folder;

        % Create MD5 Checksum file name (without extension)
        [~, base_name, ~] = fileparts(file_name);
        MD5_Checksum_file_name = strcat('MD5_', base_name, '.txt');

        checksum_file_index = find(arrayfun(@(x) strcmp(x.name, MD5_Checksum_file_name), Input_Checksum_list));
        
        if ~isempty(checksum_file_index)
            MD5_Checksum_file_path = Input_Checksum_list(checksum_file_index).folder;
            
            % Construct full path to the MD5 checksum file
            MD5_Checksum_file_path = fullfile(MD5_Checksum_file_path, MD5_Checksum_file_name);
        
            % Open and read the MD5 checksum file
            fid = fopen(MD5_Checksum_file_path, 'r');
            if fid ~= -1
                MD5_Checksum_file = fread(fid, '*char')';
                fclose(fid);
            end
        end

        % Initialize feature struct for this file (Dont'change)
        feature = struct();

        % Add file info (Dont'change)
        feature.name = file_name;
        feature.rel_path = file_path;
        feature.md5chsum = MD5_Checksum_file;
        feature.sampling_freq = Params.Fs;
        feature.mains_freq = Params.Powerline_frequency;
        feature.num_ch = Params.Number_Channels;

        %% Necessary Function (Dont'change)
        % Load the data from the file
        Num_Fun = 1;
        method = list_functions.Fun_0;
        feature.methods.(['method_' num2str(Num_Fun)]).name = method.name;
        feature.methods.(['method_' num2str(Num_Fun)]).rel_path = method.rel_path;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_md5chsum = method.codebase_md5chsum;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_git_repo = method.codebase_git_repo;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_git_commit_id = method.codebase_git_commit_id;

        % Call the function to process the data
        [result, data] = functioncaller_0(method.name, file_name, file_path, Params);

        % Store the input info of the function
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.name = 'whole_data';
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.type = 'double';
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.channels = 0;

        % Store the params of the function
        % feature.methods.(['method_' num2str(N_Fun)]).params.name = [];
        % feature.methods.(['method_' num2str(N_Fun)]).params.type = [];
        % feature.methods.(['method_' num2str(N_Fun)]).params.val = [];
                
        % Store the result of the function
        % feature.methods.(['method_' num2str(N_Fun)]).outputs.name = [];
        % feature.methods.(['method_' num2str(N_Fun)]).outputs.type = [];
        % feature.methods.(['method_' num2str(N_Fun)]).outputs.val = [];

        feature.methods.(['method_' num2str(Num_Fun)]).warning = result.warning;
        feature.methods.(['method_' num2str(Num_Fun)]).success = result.success;

        %% Part 1: Function_Num1
        Num_Fun = Num_Fun + 1;
        method = list_functions.Fun_1;
        feature.methods.(['method_' num2str(Num_Fun)]).name = method.name;
        feature.methods.(['method_' num2str(Num_Fun)]).rel_path = method.rel_path;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_md5chsum = method.codebase_md5chsum;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_git_repo = method.codebase_git_repo;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_git_commit_id = method.codebase_git_commit_id;

        channel_name = cell(1, Params.Number_Channels);
        warning = cell(1, Params.Number_Channels);
        success = zeros(1, Params.Number_Channels);
        peaks = cell(1, Params.Number_Channels);
        peak_indexes_consensus = cell(1, Params.Number_Channels);

        % Process each channel individually
        for k = 1:Params.Number_Channels
            % Call the second function to detect peaks
            result_1 = functioncaller_1(method.name, data, k, file_name, Params);
            channel_name(k) = {k};
            success(k) = result_1.success;
            warning{k} = result_1.warning;
            peaks{k} = result_1.peaks;            
            peak_indexes_consensus{k} = result_1.peak_indexes_consensus;
        end
        
        % Store the input info of the function
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.name = 'channel_data';
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.type = 'double';
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.channels = channel_name;

        % Store the params of the function
        feature.methods.(['method_' num2str(Num_Fun)]).params.name = 'fs';
        feature.methods.(['method_' num2str(Num_Fun)]).params.type = 'float';
        feature.methods.(['method_' num2str(Num_Fun)]).params.val = Params.Fs;

        % Store the result of the function
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.output_1.name = 'peaks';
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.output_1.type = 'double';
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.output_1.val = peaks;
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.output_2.name = 'peak_indexes_consensus';
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.output_2.type = 'double';
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.output_2.val = peak_indexes_consensus;
        feature.methods.(['method_' num2str(Num_Fun)]).warning = warning;
        feature.methods.(['method_' num2str(Num_Fun)]).success = success; 

        %% Part 2: Function_Num2
        Num_Fun = Num_Fun + 1;
        method = list_functions.Fun_2;
        feature.methods.(['method_' num2str(Num_Fun)]).name = method.name;
        feature.methods.(['method_' num2str(Num_Fun)]).rel_path = method.rel_path;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_md5chsum = method.codebase_md5chsum;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_git_repo = method.codebase_git_repo;
        feature.methods.(['method_' num2str(Num_Fun)]).codebase_git_commit_id = method.codebase_git_commit_id;

        % Call the third function to calculate max values
        result_2 = functioncaller_2(method.name, data, file_name, Params);

        % Store the input info of the function
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.name = 'Whole_data';
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.type = 'double';
        feature.methods.(['method_' num2str(Num_Fun)]).inputs.channels = 0;

        % Store the params of the function
        % feature.methods.(['method_' num2str(N_Fun)]).params.name = [];
        % feature.methods.(['method_' num2str(N_Fun)]).params.type = [];
        % feature.methods.(['method_' num2str(N_Fun)]).params.val = [];
                
        % Store the result of the function
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.name = 'maxvalues';
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.type = 'double';
        feature.methods.(['method_' num2str(Num_Fun)]).outputs.val = result_2.maxvalues;
        feature.methods.(['method_' num2str(Num_Fun)]).warning = result_2.warning;
        feature.methods.(['method_' num2str(Num_Fun)]).success = result_2.success;

        %% Part 3: Function_Num3
        % Your feature extarction fucntion


        %% Store the feature data in the analysis struct (Don't Change)
        analysis.(sprintf('record_%d', i)) = feature;
    
    end

    %% Save (Don't Change)

    % Convert the extracted analysis to JSON format
    json_data = jsonencode(analysis, PrettyPrint=true);
    
    % Create full path for the output JSON file
    output_file_path = fullfile(Data_Output_Path, Desired_Output_File);
    
    % Write JSON data to file
    fid = fopen(output_file_path, 'w');
    if fid == -1
        error('Cannot create JSON file %s', output_file_path);
    end
    fwrite(fid, json_data, 'char');
    fclose(fid);

    % Construct the full path to the directory
    full_path = fullfile(Data_Output_Path, "MD5_Checksum_records_matlab");

    % Check if the directory exists
    if isfolder(full_path)
        % Remove the directory and all its contents
        rmdir(full_path, 's');
    end
end