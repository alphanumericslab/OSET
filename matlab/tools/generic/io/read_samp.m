function [data, fs] = read_samp(record)
% read WFDB .dat file without WFDB library support.
% see https://physionet.org/physiotools/wag/signal-5.htm for more reference
%
% input:
%           record: WFDB data record name
% 
% output:
%           data: waveform data
%           fs: sampling frequency
% 
% Written by Qiao Li, Feb. 2025
% 
% Update by Reza Sameni, Jul. 2025 to export sampling frequency (fs)
% 
% Example - Read a record named 100.dat in current folder. 
% Note: must come with 100.hea together
% data = read_samp('100');

data=[];

% read sigInfo from the header file
sigInfo = wfdb_desc(record);


channels=sigInfo(1).channels;
data_length=sigInfo(1).LengthSamples;
sig_format=sigInfo(1).Format;

switch sig_format
    case '8'
        data_format = 'int8';
    case '16'
        data_format = 'int16';
    case '80'
        data_format = 'uint8';
    case '212'
        data_format = 'uint8';
    case '24'
        data_format = 'uint8';
    case '32'
        data_format = 'int32';
        
    otherwise
        data_format = '';
        disp('Data format currently unsupported. Exit.')
        return 
end

% read samp from .dat file
record_dat=[record '.dat'];
fileID=fopen(record_dat);
switch sig_format
    case '212'
        data=fread(fileID,[1,round(channels*data_length*1.5)],data_format);
    case '24'
        data=fread(fileID,[1,round(channels*data_length*3)],data_format);       
    otherwise
        data=fread(fileID,[channels,data_length],data_format);
end
fclose(fileID);

% process different formats
switch sig_format
    case '8'
        for i=1:size(data,1)
            raw_data(i,:) = decode_8bit_first_difference(data(i,:),sigInfo(i).InitialValue);
        end
        data=raw_data;
    case '80'
        data=data-128;
    case '212'
        data=uint8(data);
        raw_data = decode_12bit_format212(data);
        data=double(reshape(raw_data,[channels,data_length]));        
    case '24'
        data=uint8(data);
        raw_data = decode_int24_lsb(data);
        data=double(reshape(raw_data,[channels,data_length]));        
end

% process edge cases of the data
switch sig_format
%     case '80'
%         data(find(data==-128))=NaN; % ???
    case '16'
        data(data==-power(2,16-1))=NaN;
    case '24'
        data(data==-power(2,23-1))=NaN;
    case '32'
        data(data==-power(2,32-1))=NaN;
end

% process baseline and gain to restore the physical unit
sum_baseline=sum([sigInfo.Baseline]);
for i=1:channels
    if sum_baseline~=0
        data(i,:)=data(i,:)-sigInfo(i).Baseline; % -sigInfo(i).AdcZero;
    else
        data(i,:)=data(i,:)-sigInfo(i).AdcZero;
    end
    gain = sigInfo(i).Gain;
    data(i,:)=data(i,:)./gain;
end

data=data';
fs = sigInfo.SamplingFrequency;

end

function raw_data = decode_8bit_first_difference(encoded_data, initial_value)
    % Function to reconstruct raw data from 8-bit first difference encoding
    % 
    % Inputs:
    %   encoded_data - Array of 8-bit first difference values
    %   initial_value - Initial starting value from the header file
    % 
    % Output:
    %   raw_data - Reconstructed raw signal
    
    % Preallocate output array for efficiency
    raw_data = zeros(size(encoded_data));
    raw_data(1) = initial_value + encoded_data(1);
    
    % Reconstruct the raw signal by summing up first differences
    for i = 2:length(encoded_data)
        raw_data(i) = raw_data(i-1) + encoded_data(i);
    end
end

function raw_data = decode_12bit_format212(encoded_bytes)
    % Function to reconstruct raw data from 12-bit two's complement format 212
    % 
    % Input:
    %   encoded_bytes - Array of 8-bit bytes representing 12-bit samples
    % 
    % Output:
    %   raw_data - Reconstructed 16-bit vector of samples
    
    num_samples = floor(length(encoded_bytes) * 2 / 3); % Two samples per 3 bytes
    raw_data = zeros(num_samples, 1, 'int16');
    
    index = 1;
    for i = 1:3:length(encoded_bytes)-2
        % First 12-bit sample
        raw_data(index) = bitshift(int16(encoded_bytes(i+1)), 8) + int16(encoded_bytes(i));
        raw_data(index) = bitand(raw_data(index), 4095); % Mask to 12 bits
        if raw_data(index) > 2047 % Convert to signed 12-bit
            raw_data(index) = raw_data(index) - 4096;
        end
        index = index + 1;
        
        % Second 12-bit sample
        raw_data(index) = bitshift(int16(bitand(uint8(encoded_bytes(i+1)),0xF0)), 4) + int16(encoded_bytes(i+2));
        raw_data(index) = bitand(raw_data(index), 4095); % Mask to 12 bits
        if raw_data(index) > 2047 % Convert to signed 12-bit
            raw_data(index) = raw_data(index) - 4096;
        end
        index = index + 1;
    end
end

function int24_values = decode_int24_lsb(byte_data)
    % Convert 3-byte sequences into signed 24-bit integers (little-endian)
    % byte_data: Nx3 matrix, where each row is a 24-bit integer in LSB-first format

    % Ensure byte_data is a 3-column matrix
    if mod(numel(byte_data), 3) ~= 0
        error('Input data length must be a multiple of 3.');
    end

    % Reshape into Nx3 (each row is a 24-bit integer)
    byte_data = reshape(byte_data, 3, [])';

    % Convert to int32 (to avoid overflow)
    int24_values = int32(byte_data(:,1)) + ...
                   bitshift(int32(byte_data(:,2)), 8) + ...
                   bitshift(int32(byte_data(:,3)), 16);

    % Handle negative values (two’s complement sign extension)
    negative_mask = bitshift(int32(1), 23); % 0x00800000
    int24_values(int24_values >= negative_mask) = int24_values(int24_values >= negative_mask) - bitshift(int32(1), 24); % -2^24
end