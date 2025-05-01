function [data, fs, sigInfo] = rdmat_wfdb(varargin)
% Description: Import WFDB signals from .mat files in physical units
%
% This function is part of the OSET (Open-Source Electrophysiological Toolbox)
% and reads WFDB format signals from .mat files that were previously
% converted using WFDB2MAT. It handles signal conversion from digital to physical
% units and provides signal information from the associated header file.
%
% INPUT:
% recordName - Name of the WFDB record (without extension)
%              Format: string
%              Example: '200m' (for files 200m.mat and 200m.hea)
%              Note: Both .mat and .hea files must be present
%
% OUTPUT:
% data - Signal data in physical units
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
%           - Units: physical units
%           - Baseline: baseline value
%           - Gain: signal gain
%           - Description: channel description
%           - fmt: WFDB format code
%
% PROCESSING STEPS:
% 1. Read header file (.hea) to get signal information
% 2. Load signal data from .mat file
% 3. Convert digital values to physical units:
%    - Handle format-specific digital values (80, 16, etc.)
%    - Apply baseline correction
%    - Apply gain scaling
% 4. Generate time vector
%
% NOTE:
% - This function is significantly faster than RDSAMP as it reads pre-converted
%   .mat files instead of processing raw WFDB files
% - For integrity checking, compare output with RDSAMP:
%   [data1, fs1, info1] = rdmat_wfdb('200m');
%   [data2, fs2] = rdsamp('200m');
%   if sum(abs(data1-data2)) ~= 0
%       error('Record not compatible with rdmat_wfdb');
%   end
%
% LIMITATIONS:
% - Does not support all WFDB record format features
% - Does not support multiresolution signals
%
% Author: Ikaro Silva
% Date: November 26, 2014
% Version: 1.2
%
% Example:
%   wfdb2mat('mitdb/200');
%   [data, fs, info] = rdmat_wfdb('200m');
%
% See also RDSAMP, WFDB2MAT

%endOfHelp

%Set default pararameter values
inputs={'recordName'};
defGain=200; %Default value for missing gains

for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};']);
    end
end


fid = fopen([recordName, '.hea'], 'rt');
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

load([recordName '.mat']);

for m = 1:M
    % Interpreting digital values of byte offset format 80
    if strcmp(siginfo(m).fmt, '80') 
        val(m,:)=val(m,:)-128;
        wfdbNaN=-128;
    elseif strcmp(siginfo(m).fmt, '16')
        wfdbNaN=-32768;
    else
        wfdbNaN=-2147483648;
    end
    
    % Fill in NaNs before subtracting and dividing. 
    val(m, val(m,:)==wfdbNaN)=nan;
    
    %Convert from digital units to physical units.
    % Mapping should be similar to that of rdsamp.c:
    % http://www.physionet.org/physiotools/wfdb/app/rdsamp.c
    val(m, :) = (val(m, :) - siginfo(m).Baseline ) / siginfo(m).Gain;
end

%Reshape to the Toolbox's standard format
val=val'; 

%Generate time vector
N=size(val,1);
tm =linspace(0,(N-1)/Fs,N);

data = val;
fs = Fs;
sigInfo = siginfo;

end
