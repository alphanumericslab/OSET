function feature_handle_csv(output_csv_file_name, features_matrix, feature_names, features_units, feature_description,  fs, window_time_info, source_names, num_decimal_places, overwrite_flag)

% Description: write extracted features in output_csv_filename file with standard OSET AI-ML ready format
%
% INPUT:
% output_csv_filename -
% features_matrix   -
% feature_names -
% features_units -
% feature_description -
% fs -
% window_time_info -
% num_decimal_places -
% overwrite_flag -

% Author: Sajjad Karimi
% Date: Mar 20, 2025
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com

if nargin<4 || isempty(features_units)
    feature_names = repmat({'NA'}, 1, size(features_matrix,2));
end

if nargin<4 || isempty(features_units)
    features_units = repmat({'scaler'}, 1, size(features_matrix,2));
end

if nargin<5 || isempty(feature_description)
    feature_description = repmat({"NA"}, 1, size(features_matrix,2));
end

if nargin<6 || isempty(fs)
    fs = 250;
end

if nargin<7 || isempty(window_time_info)
    window_time_info = [1,fs*10];
end

if nargin<8 || isempty(source_names)
    source_names = repmat({"NA"}, size(features_matrix,1),1);
end

if nargin<9 || isempty(num_decimal_places)
    num_decimal_places = 6;
end

if nargin<10 || isempty(overwrite_flag)
    overwrite_flag = false;
end


if overwrite_flag && exist(output_csv_file_name, 'file')
    delete(output_csv_file_name); % Delete the file
    fprintf('File "%s" deleted.\n', output_csv_file_name);
end

% Check if file does not exist make .csv file with headers, units, gains
% and description
if ~exist(output_csv_file_name, 'file')

    % Get the folder path from the file name
    output_folder = fileparts(output_csv_file_name);

    % Check if the folder exists, and if not, create it
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    feature_names = feature_names(:).';  % Convert to row vector;
    features_units = features_units(:).';
    feature_description = feature_description(:).';

    % Header information
    csv_var_names = [{'source', 'start_time', 'stop_time'}, feature_names];

    csv_var_description = [{"Indicate source channel or modality for each row", "Starting index of this window", "Stoping index of this window"}, feature_description];

    csv_units = [{'na','s','s'}, features_units ];

    csv_gain = [{'1',num2str(1/fs),num2str(1/fs)}, repmat({'1'}, 1, length(feature_names))];


    % Open the file for writing
    fid = fopen(output_csv_file_name, 'w');

    % Write the variable names
    fprintf(fid, '%s,', csv_var_names{1:end-1});
    fprintf(fid, '%s\n', csv_var_names{end});

    % Write the header descriptions
    fprintf(fid, '%s,', csv_var_description{1:end-1});
    fprintf(fid, '%s\n', csv_var_description{end});

    % Write the units
    fprintf(fid, '%s,', csv_units{1:end-1});
    fprintf(fid, '%s\n', csv_units{end});

    % Write the gain values
    fprintf(fid, '%s,', csv_gain{1:end-1});
    fprintf(fid, '%s\n', csv_gain{end});

    % Close the file
    fclose(fid);

end

first_index = window_time_info(:,1);
stop_index = window_time_info(:,2);

csv_data = [first_index, stop_index, features_matrix];

% Limit the number of decimal places to reduce the file sizes
csv_data = round(csv_data * 10^num_decimal_places) / 10^num_decimal_places;

% Append the data to the file
% writematrix(csv_data, output_csv_file_name, 'WriteMode', 'append');

csv_data = [source_names, num2cell(csv_data)];
writecell(csv_data, output_csv_file_name, 'WriteMode', 'append');


