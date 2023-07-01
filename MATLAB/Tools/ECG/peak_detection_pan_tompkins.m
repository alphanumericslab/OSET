function [peaks, peak_indexes] = peak_detection_pan_tompkins(data, fs, varargin)
% peak_detection_pan_tompkins - R-peak detector based on Pan-Tompkins method.
%
%   [peaks, peak_indexes] = peak_detection_pan_tompkins(data, fs, varargin)
%
%   This function implements the Pan-Tompkins algorithm for R-peak detection
%   in ECG signals.
%
%   Inputs:
%       data: Vector of input ECG data
%       fs: Sampling rate in Hz
%       varargin: Optional arguments (variable number and order)
%           - varargin{1}: Struct containing algorithm parameters (optional)
%           - varargin{2:end}: Individual parameter-value pairs (optional)
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Reference:
%       Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
%       Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532
%
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Set default parameter values
default_params.fc_lp = 0.5;
default_params.fc_hp = 5;
default_params.window_length = 0.150;
default_params.threshold_ratio = 0.2;
default_params.refractory_period = 0.200;

% Parse input arguments
params = default_params;
if ~isempty(varargin) && isstruct(varargin{1})
    % params struct is provided
    params = fill_default_params(varargin{1}, default_params);
end


% Low-pass filter
[b_lp, a_lp] = butter(5, params.fc_lp / (fs/2), 'low');
filtered_data_lp = filtfilt(b_lp, a_lp, data);

% High-pass filter
[b_hp, a_hp] = butter(5, params.fc_hp / (fs/2), 'high');
filtered_data_hp = filtfilt(b_hp, a_hp, filtered_data_lp);

% Differentiation to enhance R-peaks
diff_data = diff(filtered_data_hp);

% Squaring to further emphasize R-peaks
squared_data = diff_data .^ 2;

% Moving average integration
window_length = round(params.window_length * fs);
window = ones(1, window_length) / window_length;
integrated_data = conv(squared_data, window, 'same');

% Find R-peaks using adaptive thresholding
threshold = params.threshold_ratio * max(integrated_data);
peaks = integrated_data > threshold;

% Refractory period to avoid detecting multiple peaks within a short duration
refractory_period = round(params.refractory_period * fs);
peaks = remove_nearby_peaks(peaks, refractory_period);

% Get peak indexes
peak_indexes = find(peaks);

end

function peaks = remove_nearby_peaks(peaks, refractory_period)
% Helper function to remove nearby peaks within the refractory period

last_peak = find(peaks, 1);
for i = (last_peak + 1):length(peaks)
    if peaks(i)
        if i - last_peak <= refractory_period
            peaks(last_peak:i) = false;
        else
            last_peak = i;
        end
    end
end

end

function params = fill_default_params(params, default_params)
% Helper function to fill missing fields in params struct with default values

fields = fieldnames(default_params);
for i = 1:numel(fields)
    if ~isfield(params, fields{i})
        params.(fields{i}) = default_params.(fields{i});
    end
end

end
