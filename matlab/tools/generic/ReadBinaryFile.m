function [x, fs, status] = ReadBinaryFile(ifname)
%
% [x, fs, status] = ReadBinaryFile(ifname)
%
% output:
% ifname: input file name (with path and extension)
%
% outputs:
% x: matrix of input data (channels x samples)
% fs: sampling frequency
% status: 0 (write successful), -2 (file not opened), -1 (file size doesn not match header)
%
% Binary file read format:
%   fs (8 bytes), 'double', 'ieee-le'
%   rows (8 bytes), 'uint64', 'ieee-le'
%   cols (8 bytes), 'uint64', 'ieee-le'
%   data (rows x cols )x 8 bytes, 'double', 'ieee-le' (written column wise)
%
% Open Source ECG Toolbox, version 3.14, January 2020
% Copyright Reza Sameni
% reza.sameni@gmail.com

fs = 0;
x = [];
fileID = fopen(ifname, 'r');
if(fileID >= 0)
    fs = fread(fileID, 1, 'double', 'ieee-le'); % 8 bytes
    rows = fread(fileID, 1, 'uint64', 'ieee-le'); % 8 bytes
    cols = fread(fileID, 1, 'double', 'ieee-le'); % 8 bytes
    data = fread(fileID, 'double', 'ieee-le'); % (rows x cols )x 8 bytes, written column wise
    fclose(fileID);
    if( numel(data) == rows * cols )
        x = reshape(data, rows, cols);
        status = 0;
    else
        status = -1;
    end
else
    status = -2;
end
