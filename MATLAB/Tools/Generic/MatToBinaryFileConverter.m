function [status, count] = MatToBinaryFileConverter(x, fs, ofname)
%
% [status, count] = MatToBinaryFileConverter(x, fs, ofname)
%
% inputs:
% x: matrix of input data (channels x samples)
% fs: sampling frequency
% ofname: output file name (with path and extension)
%
% output:
% status: 0 (write successful), -2 (file not opened), -1 (file not closed)
% count: number of 8-byte elements written in the file
% 
% Binary file write format:
%   fs (8 bytes), 'double', 'ieee-le'
%   rows (8 bytes), 'uint64', 'ieee-le'
%   cols (8 bytes), 'uint64', 'ieee-le'
%   data (rows x cols )x 8 bytes, 'double', 'ieee-le' (written column wise)
%
% Open Source ECG Toolbox, version 3.14, January 2020
% Copyright Reza Sameni
% reza.sameni@gmail.com

[rows, cols] = size(x);
fileID = fopen(ofname, 'w+');
if(fileID >= 0)
    count = fwrite(fileID, fs, 'double', 'ieee-le'); % 8 bytes
    count = count + fwrite(fileID, rows, 'uint64', 'ieee-le'); % 8 bytes
    count = count + fwrite(fileID, cols, 'double', 'ieee-le'); % 8 bytes
    count = count + fwrite(fileID, x, 'double', 'ieee-le'); % (rows x cols )x 8 bytes, written column wise
    status = fclose(fileID);
else
    count = 0;
    status = -2;
end
