function [data, fs, sigInfo] = rddat_wfdb(record)
% Description: Read WFDB .dat file without WFDB library support
%
% This function is part of the OSET (Open-Source Electrophysiological Toolbox)
% and reads and decodes WFDB format signal files (.dat) along with their
% associated header files (.hea). It supports various WFDB formats including 8-bit,
% 16-bit, 24-bit, and format 212 (12-bit packed).
%
% INPUT:
% record - WFDB record name (without extension)
%          Format: string
%          Example: '100' (for files 100.dat and 100.hea)
%          Note: The .hea file must be present in the same directory
%
% OUTPUT:
% data - Decoded signal data
%        Format: double matrix [T x C]
%        T: number of time samples
%        C: number of channels
%        Units: physical units as specified in the header file
%
% fs - Sampling frequency in Hz
%
% sigInfo - Signal information structure
%           Format: struct array [1 x C]
%           Fields:
%           - RecordName: name of the record
%           - channels: number of channels
%           - SamplingFrequency: sampling rate in Hz
%           - LengthSamples: number of samples
%           - Format: WFDB format code
%           - Gain: signal gain
%           - Baseline: baseline value
%           - Units: physical units
%           - Description: channel description
%
% SUPPORTED FORMATS:
% '8'  - 8-bit first difference
% '16' - 16-bit integer
% '80' - 8-bit offset binary
% '212'- 12-bit format 212 (packed)
% '24' - 24-bit integer
% '32' - 32-bit integer
%
% PROCESSING STEPS:
% 1. Read header file (.hea) to get signal information
% 2. Determine data format and allocate appropriate data type
% 3. Read binary data from .dat file
% 4. Decode data based on format:
%    - 8-bit: first difference decoding
%    - 12-bit: format 212 unpacking
%    - 24-bit: little-endian conversion
% 5. Handle edge cases (NaN values)
% 6. Apply baseline correction and gain
%
% DEPENDENCIES:
% - wfdb_desc: reads and parses WFDB header files
% - decode_8bit_first_difference: decodes 8-bit first difference format
% - decode_12bit_format212: decodes 12-bit format 212
% - decode_int24_lsb: decodes 24-bit little-endian format
%
% Author: Qiao Li
% Date: Feb. 2025
% Modified by: Sajjad Karimi
% Date: Apr. 2025
%
% Example:
%   [data, info] = rddat_wfdb('100');
%   % Reads signal from 100.dat and header from 100.hea

data=[];

% read sigInfo from the header file
sigInfo = wfdb_desc(record);

% sigInfo = wfdb_info(record);

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
        raw_data = zeros(size(data,1),size(data,2));
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
fs=sigInfo(1).SamplingFrequency;

end

function raw_data = decode_8bit_first_difference(encoded_data, initial_value)
% Description: Reconstruct raw data from 8-bit first difference encoding
%
% This function decodes data that has been encoded using first difference
% encoding, where each sample represents the difference from the previous sample.
%
% INPUT:
% encoded_data - Array of 8-bit first difference values
%                Format: int8 array
%
% initial_value - Starting value for reconstruction
%                 Format: double scalar
%
% OUTPUT:
% raw_data - Reconstructed signal
%            Format: double array
%            Same length as encoded_data

% Preallocate output array for efficiency
raw_data = zeros(size(encoded_data));
raw_data(1) = initial_value + encoded_data(1);

% Reconstruct the raw signal by summing up first differences
for i = 2:length(encoded_data)
    raw_data(i) = raw_data(i-1) + encoded_data(i);
end
end

function raw_data = decode_12bit_format212(encoded_bytes)
% Description: Reconstruct raw data from 12-bit format 212
%
% This function decodes the WFDB format 212, which packs two 12-bit samples
% into three bytes. The format is commonly used in ECG recordings.
%
% INPUT:
% encoded_bytes - Array of 8-bit bytes representing 12-bit samples
%                 Format: uint8 array
%                 Length must be multiple of 3
%
% OUTPUT:
% raw_data - Reconstructed 12-bit samples
%            Format: int16 array
%            Length: floor(length(encoded_bytes) * 2/3)

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
% Description: Convert 3-byte sequences into signed 24-bit integers
%
% This function decodes 24-bit integers stored in little-endian format
% (least significant byte first).
%
% INPUT:
% byte_data - Array of bytes representing 24-bit integers
%             Format: uint8 array
%             Length must be multiple of 3
%
% OUTPUT:
% int24_values - Decoded 24-bit integers
%                Format: int32 array
%                Length: length(byte_data)/3

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

% Handle negative values (two's complement sign extension)
negative_mask = bitshift(int32(1), 23); % 0x00800000
int24_values(int24_values >= negative_mask) = int24_values(int24_values >= negative_mask) - bitshift(int32(1), 24); % -2^24
end


function sigInfo = wfdb_desc(headname)
% Description: Read and parse WFDB header file
%
% This function reads a WFDB header file (.hea) and extracts signal
% information including format, sampling rate, gain, and channel descriptions.
%
% INPUT:
% headname - Header file name without extension
%            Format: string
%
% OUTPUT:
% sigInfo - Signal information structure
%           Format: struct array
%           Fields:
%           - RecordName: record name
%           - channels: number of channels
%           - SamplingFrequency: sampling rate in Hz
%           - LengthSamples: number of samples
%           - Format: data format code
%           - Gain: signal gain
%           - Baseline: baseline value
%           - Units: physical units
%           - Description: channel description

sigInfo = struct;

fid = fopen([headname '.hea']);
str = fgets(fid);
cc = textscan(str,'%s');
cc = cc{1};

var_name = {'RecordName','channels','SamplingFrequency','LengthSamples','time', 'date'};
var_type = {'s','d','d','d','s','s'};
n=0;
for i=1:length(cc)
    if i<=length(var_name)
        com = sprintf('%s = "%s";', var_name{i}, cc{i});
        eval(com);
        n=n+1;
    end
end

item_name = {'File','Format','Gains','AdcResolution','AdcZero','InitialValue','CheckSum','BlockSize','Description'};
item_type = {'s','s','s','d','d','d','d','d','s'};
n_files=0;
n_item=0;
str = fgets(fid);

[~, filename, ~] = fileparts(headname);
while ischar(str)

    if length(str)>length(filename)

        if  ~strcmp(str(1:length(filename)),filename)
            str = fgets(fid);
            continue
        end
    else
        str = fgets(fid);
        continue
    end

    cc2 = textscan(str,'%s');
    cc2 = cc2{1};
    n_files=n_files+1;
    for i=1:n
        com = sprintf('sigInfo(%d).%s = "%s";',n_files, var_name{i}, cc{i});
        eval(com);
        if var_type{i} == 'd'
            com = sprintf(['sigInfo(%d).%s = str2num(sigInfo(%d).%s);'],n_files,var_name{i},n_files,var_name{i});
            eval(com);
        end
    end
    if n>=6
        sigInfo(n_files).StartTime=strjoin(['[', sigInfo(n_files).time, ' ', sigInfo(n_files).date, ']']);
    else
        sigInfo(n_files).StartTime="";
    end
    for i=1:length(cc2)
        if i<=length(item_name)
            com = sprintf(['sigInfo(%d).%s = "%s";'],n_files,item_name{i},cc2{i});
            eval(com);
            if item_type{i} == 'd'
                com = sprintf(['sigInfo(%d).%s = str2num(sigInfo(%d).%s);'],n_files,item_name{i},n_files,item_name{i});
                eval(com);
            end
            n_item = n_item + 1;
        end
    end

    numStr = regexp(sigInfo(n_files).Gains, '[+-]?\d+(\.\d+)?', 'match');

    sigInfo(n_files).Gain=str2num(numStr{1});

    if length(numStr)>1
        sigInfo(n_files).Baseline=str2num(numStr{2});
    else
        sigInfo(n_files).Baseline = 0;
    end

    k=strfind(sigInfo(n_files).Gains, '/');
    sigInfo(n_files).Units=sigInfo(n_files).Gains{1}(k+1:end);

    str = fgets(fid);
end
fclose(fid);
end



function siginfo = wfdb_info(headname)

fid = fopen([headname, '.hea'], 'rt');
if(fid==-1)
    error(['Could not open file: ' recordName '.hea !']);
end

%Following the documentation described in :
%http://www.physionet.org/physiotools/wag/header-5.htm
%to parse the header file

%Skip any comment lines
str=fgetl(fid);
while(strcmp(str(1),'#'))
    str=fgetl(fid);
end

%Process Record Line Info
info=textscan(str,'%s %d %f %d %s %s');
M=info{2}; %Number of signals present
Fs=info{3};

%Process Signal Specification lines. Assumes no comments between lines.
siginfo(M)=struct();
for m = 1:M
    str=fgetl(fid);
    info=textscan(str,'%s %s %s %d %d %f %d %d %s');
    fmt=info{2}{:};
    gain=info{3}{:};

    %Get Signal Units if present
    ind=strfind(gain,'/');
    if(~isempty(ind))
        siginfo(m).Units=gain(ind+1:end);
        gain=gain(1:ind-1);
    end

    %Get Signal Baseline if present
    ind=strfind(gain,'(');
    if(~isempty(ind))
        ind2=strfind(gain,')');
        siginfo(m).Baseline=str2num(gain(ind+1:ind2-1));
        gain=gain(1:ind-1);
    else
        %If Baseline is missing, set it equal to ADC Zero
        adc_zero=info{5};
        if(~isempty(adc_zero))
            siginfo(m).Baseline=double(adc_zero);
        else
            error('Could not obtain signal baseline');
        end
    end

    %Get Signal Gain
    gain=str2num(gain);
    if(gain==0)
        %Set gain to default value in this case
        gain=defGain;
    end
    siginfo(m).Gain=double(gain);

    %Get Signal Descriptor
    siginfo(m).Description=info{9}{:};

    % Store format for later
    siginfo(m).fmt=fmt(1:strfind(fmt,'+')-1);

end
fclose(fid);

end