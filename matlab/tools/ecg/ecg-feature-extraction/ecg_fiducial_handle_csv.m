function ecg_fiducial_handle_csv(output_csv_file_name, ecg_data, fs, ecg_fiducial_position, lead_names, first_index, overwrite_flag)

% all_features_vector = extract_ecg_features_one_record(signal_name, data, fs, n_eigenvalues, n_features, f_notch)
% Description: Extract features from one record of ECG signal
%
% INPUT:
% ecg_data - ECG signal as matrix CxT which C indicates number of ECG channels (in microvolts).
% fs   - Sampling frequency of the ECG signal (in Hz).
% lead_names - a cell array containing each ECG lead name
% n_svd - Number of eigenvalues for SVD features (default 5).
% num_morpho_samples - Number of samples for ECG average beat (default 25).
% norm_flag - A boolian value for normalizing ECG average beat before sampling (default true).
%
% DEPENDENCIES:
% 1. peak_det_likelihood function from the OSET package
% 2. peak_det_likelihood_long_recs function from the OSET package
% 3. fiducial_det_lsim function from the OSET package
%
% OUTPUT:
% ecg_features_vector - includes the following features per each channel:
%   1. SNR features (2)
% vecg_feature_info - includes following information feilds
%   ecg_feature_names - name of features for all channels
%   ecg_features_units - unit of all extracted features
%   ecg_fiducial_position - cell array of detected fiducial points of ECG
% exit_flag - with 0 indicates successful, and minus values (-K)
%   indicates error in processing K channels out of C channels

% Author: Sajjad Karimi
% Date: Mar 20, 2025
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com


if nargin<7 || isempty(overwrite_flag)
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

    feature_names = {'p_on', 'p', 'p_off', 'qrs_on', 'rpeak_indexes', 'qrs_off', 't_on', 't', 't_off',...
        'p_on_amp', 'p_amp', 'p_off_amp', 'qrs_on_amp', 'r_amp', 'qrs_off_amp', 't_on_amp', 't_amp', 't_off_amp' };

    feature_description = {"P-wave onset index", "P-wave peak index", "P-wave offset index", ...
        "QRS onset index", "R-peak index (from delineator)", "QRS offset index", ...
        "T-wave onset index", "T-wave peak index", "T-wave offset index", ...
        "P-wave onset amplitude", "P-wave peak amplitude", "P-wave offset amplitude", ...
        "QRS onset amplitude", "R-peak amplitude", "QRS offset amplitude", ...
        "T-wave onset amplitude", "T-wave peak amplitude", "T-wave offset amplitude"};


    features_units =  [repmat({'s'},1,length(feature_names)/2),  repmat({'mv'},1,length(feature_names)/2)];


    % % Header information
    % csv_var_names = [{'lead_name'}, feature_names];
    % csv_var_description = [{"Indicate which lead fiducial points are reported"}, feature_description];

    % Concatenate with commas
    stacked_leads = strjoin(lead_names, ' ');
    % Header information
    csv_var_names = [{'lead_num'}, feature_names];
    % csv_var_description = [{strjoin([string(stacked_leads), "(Indicate lead number for reported fiducial points)"])}, feature_description];
    csv_var_description = [{string(stacked_leads)}, feature_description];

    csv_units = [{'NA'}, features_units ];

    csv_gain = [{'1'}, repmat({num2str(1/fs)},1,length(feature_names)/2),  repmat({'1'},1,length(feature_names)/2)];


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


C = length(ecg_fiducial_position);
num_decimal_places = 6; % for Amplititude values 

for c = 1:C

    fiducials = ecg_fiducial_position{c};

    features_matrix_index =  [fiducials.Pon(:), fiducials.P(:), fiducials.Poff(:), fiducials.QRSon(:), fiducials.R(:), fiducials.QRSoff(:), fiducials.Ton(:), fiducials.T(:), fiducials.Toff(:)];
    ecg_data_lead = ecg_data(c,:);
    features_matrix_index_fill = features_matrix_index;
    features_matrix_index_fill(isnan(features_matrix_index_fill)) = 1;


    features_matrix_amp = ecg_data_lead(features_matrix_index_fill);
    % Limit the number of decimal places to reduce the file sizes
    features_matrix_amp = round(features_matrix_amp * 10^num_decimal_places) / 10^num_decimal_places;

    features_matrix_amp(isnan(features_matrix_index)) = nan;


    csv_data = [ repmat(c,size(features_matrix_index,1),1), (-1 + first_index + features_matrix_index) , features_matrix_amp];

    % Append the data to the file
    writematrix(csv_data, output_csv_file_name, 'WriteMode', 'append');


end
