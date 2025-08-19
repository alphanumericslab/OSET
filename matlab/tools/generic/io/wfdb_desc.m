function sigInfo = wfdb_desc(headname)
%WFDB_DESC  Read and parse a WFDB .hea header into a struct array.
%
%   sigInfo = wfdb_desc('recordname')
%
%   Fields produced per channel (mirroring your original intent):
%     RecordName, channels, SamplingFrequency, LengthSamples, time, date, StartTime
%     File, Format, Gains, AdcResolution, AdcZero, InitialValue, CheckSum, BlockSize, Description
%     Gain, Baseline, Units

sigInfo = struct([]);

fname = [headname '.hea'];
[fid, msg] = fopen(fname, 'r');
if fid == -1
    error('wfdb_desc:openFailed', 'Cannot open %s: %s', fname, msg);
end
c = onCleanup(@() fclose(fid));

% --- First line (record summary) ------------------------------------------
firstLine = fgetl(fid);
if ~ischar(firstLine)
    error('wfdb_desc:emptyHeader', 'Header file %s is empty.', fname);
end
cc = strsplit(strtrim(firstLine));

% Expected tokens: recordName, nsig, fs, length, [date time ...]
if numel(cc) < 3
    error('wfdb_desc:badHeader', 'Header summary line is incomplete in %s.', fname);
end

RecordName        = cc{1};
channels          = str2double(cc{2});
SamplingFrequency = str2double(cc{3});

% LengthSamples is sometimes absent or replaced by duration; guard it
if numel(cc) >= 4 && ~startsWith(cc{4}, '#')
    LengthSamples = str2double(cc{4});
else
    LengthSamples = NaN;
end

% Anything beyond 4th token we treat as date/time blob (WFDB variants exist)
timeStr = '';
dateStr = '';
if numel(cc) >= 6
    % Your original order used time then date; keep that but guard types
    timeStr = cc{5};
    dateStr = strjoin(cc(6:end), ' ');
elseif numel(cc) == 5
    % Only one token given (either date or time); store as "time"
    timeStr = cc{5};
end

% Preallocate struct array for channels if channels is valid; otherwise grow
if ~isfinite(channels) || channels < 1
    channels = 0; % we’ll grow dynamically if needed
end

% Helper to stamp common fields into a given element
function s = stamp_common(s)
    s.RecordName        = RecordName;
    s.channels          = channels;
    s.SamplingFrequency = SamplingFrequency;
    s.LengthSamples     = LengthSamples;
    s.time              = timeStr;
    s.date              = dateStr;
    if ~isempty(timeStr) && ~isempty(dateStr)
        s.StartTime = sprintf('[%s %s]', timeStr, dateStr);
    else
        s.StartTime = '';
    end
end

% --- Per-signal lines ------------------------------------------------------
sigInfo = repmat(stamp_common(struct()), max(channels,1), 1);
n_files = 0;

line = fgetl(fid);
while ischar(line) && ~(~isempty(line) && line(1) == '#')
    str = strtrim(line);
    if isempty(str)
        line = fgetl(fid);
        continue; % skip blank lines
    end

    toks = strsplit(str);
    % Map to your item fields by position (brittle but mirrors original)
    item_name = {'File','Format','Gains','AdcResolution','AdcZero','InitialValue','CheckSum','BlockSize','Description'};
    item_type = {'s','s','s','d','d','d','d','d','s'};

    n_files = n_files + 1;
    if n_files > numel(sigInfo)
        % Grow if channels was missing/wrong
        sigInfo(end+1,1) = stamp_common(struct());
    end

    % Stamp common fields
    sigInfo(n_files) = stamp_common(sigInfo(n_files));

    % Fill items by position safely
    for i = 1:min(numel(toks), numel(item_name))
        fld = item_name{i};
        val = toks{i};
        switch item_type{i}
            case 'd'
                numv = str2double(val);
                if ~isfinite(numv)
                    numv = NaN;
                end
                sigInfo(n_files).(fld) = numv;
            otherwise
                sigInfo(n_files).(fld) = val;
        end
    end

    % Ensure present fields exist even if missing in line
    for i = 1:numel(item_name)
        fld = item_name{i};
        if ~isfield(sigInfo(n_files), fld)
            if item_type{i} == 'd'
                sigInfo(n_files).(fld) = NaN;
            else
                sigInfo(n_files).(fld) = '';
            end
        end
    end

    % --- Derive Gain / Baseline / Units from Gains field -------------------
    gainsField = sigInfo(n_files).Gains; % as char
    gainVal   = NaN;
    baseVal   = NaN;
    unitsVal  = '';

    if ~isempty(gainsField)
        % Extract first two numbers: gain and baseline (if present)
        numStr = regexp(gainsField, '[+-]?\d+(\.\d+)?', 'match');
        if ~isempty(numStr)
            gainVal = str2double(numStr{1});
            if numel(numStr) >= 2
                baseVal = str2double(numStr{2});
            end
        end
        % Units after first '/' if present
        k = strfind(gainsField, '/');
        if ~isempty(k)
            unitsVal = strtrim(gainsField(k(1)+1:end));
        end
    end
    sigInfo(n_files).Gain     = gainVal;
    sigInfo(n_files).Baseline = baseVal;
    sigInfo(n_files).Units    = unitsVal;

    line = fgetl(fid);
end
end
