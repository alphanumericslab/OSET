clear
close all
clc

fname = 'ICARE_0097_20000102_185034_segment';
load(fname, 'signal', 'sampling_rate', 'channel_names');
signal = signal(:, 1:10000);

save('sample.mat', 'signal', 'sampling_rate', 'channel_names'); 

% Save the raw data matrix in CSV format
fout_csv = 'sample_raw.csv';
writematrix(signal', fout_csv);
% Read back the file
signal_text = num2cell(readmatrix(fout_csv));
% Append the channel names and units in the headers
combined_data = string([channel_names ; repmat({'mV'}, 1, size(signal, 1)) ; signal_text]);
% Write the combined data back to the CSV file
writematrix(combined_data, fout_csv);

