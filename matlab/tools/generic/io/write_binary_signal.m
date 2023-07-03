function [status, count] = write_binary_signal(x, fs, ofname)
% write_binary_signal - Save a multichannel signal as a binary file (readable by read_binary_signal)
%
% Syntax: [status, count] = write_binary_signal(x, fs, ofname)
% 
% Inputs:
%   x: Matrix of input data (channels x samples)
%   fs: Sampling frequency
%   ofname: Output file name (including path and extension)
%
% Outputs:
%   status: Status of the file write operation
%           0 - Write successful
%          -2 - File not opened
%          -1 - File not closed
%   count: Number of 8-byte elements written in the file
%
% Binary file write format:
%   fs (8 bytes), 'double', 'ieee-le'
%   rows (8 bytes), 'uint64', 'ieee-le'
%   cols (8 bytes), 'uint64', 'ieee-le'
%   data (rows x cols) x 8 bytes, 'double', 'ieee-le' (written column-wise)
%
% Revision History:
%   2020: First release
%   2023: Renamed from deprecated version MatToBinaryFileConverter()
%
% References:
%   Reza Sameni, 2020-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

[rows, cols] = size(x); % Get the dimensions of the input data
fileID = fopen(ofname, 'w+'); % Open the file for writing
if fileID >= 0
    % Write the data to the file
    count = fwrite(fileID, fs, 'double', 'ieee-le'); % Write sampling frequency (8 bytes)
    count = count + fwrite(fileID, rows, 'uint64', 'ieee-le'); % Write number of rows (8 bytes)
    count = count + fwrite(fileID, cols, 'double', 'ieee-le'); % Write number of columns (8 bytes)
    count = count + fwrite(fileID, x, 'double', 'ieee-le'); % Write data (rows x cols) x 8 bytes, written column-wise
    status = fclose(fileID); % Close the file
else
    count = 0;
    status = -2; % File not opened
end
