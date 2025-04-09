
clc
clear
close all


input_wfdb_address = './sample-ecg/HR00001.mat';
ecg_csv_file_name = './sample-ecg/HR00001_ecg.csv';

lead_names_target = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
window_dur = 10;
start_time = 1;
stop_time = 0;


process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time);

input_wfdb_address = './sample-ecg/sel30.dat';
ecg_csv_file_name = './sample-ecg/sel30_ecg.csv';

lead_names_target = {'ECG1', 'ECG2', 'V4', 'V5', 'V6'};
window_dur = 60;
start_time = 10;
stop_time = 430;


process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time);