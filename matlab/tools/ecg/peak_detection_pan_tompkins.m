function [peaks, peak_indexes] = peak_detection_pan_tompkins(data, fs, varargin)
% peak_detection_pan_tompkins - R-peak detector based on Pan-Tompkins method.
%
%   [peaks, peak_indexes] = peak_detection_pan_tompkins(data, fs, varargin)
%
%   This function implements the Pan-Tompkins algorithm for R-peak detection
%   in ECG signals, with a simplified post-detection R-peak selection logic
%
%   Inputs:
%       data: Vector of input ECG data
%       fs: Sampling rate in Hz
%       fc_low (Optional): BP filter lower cutoff frequency in Hz (default: 5.0 Hz)
%       fc_high (Optional): BP filter upper cutoff frequency in Hz (default: 15.0 Hz)
%       window_length (Optional): Integration window length in seconds (default: 0.150 s)
%       threshold_ratio (Optional): Threshold ratio for peak detection (default: 0.2)
%       refractory_period (Optional): Refractory period in seconds (default: 0.2 s)
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Reference:
%       Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
%       Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532
%
%   Reza Sameni, 2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Set default parameter values
default_params.fc_low = 5.0;
default_params.fc_high = 15;
default_params.window_length = 0.150;
default_params.threshold_ratio = 0.2;
default_params.refractory_period = 0.200;

% Parse input arguments
params = default_params;
if ~isempty(varargin) && isstruct(varargin{1})
    % params struct is provided
    params = fill_default_params(varargin{1}, default_params);
end

data = data(:);

% High-pass filter
[b_hp, a_hp] = butter(5, params.fc_low / (fs/2), 'high');
filtered_data_hp = filtfilt(b_hp, a_hp, data);

% Low-pass filter
[b_lp, a_lp] = butter(5, params.fc_high / (fs/2), 'low');
filtered_data = filtfilt(b_lp, a_lp, filtered_data_hp);

% Differentiation to enhance R-peaks
diff_data = [0 ; diff(filtered_data)];

% Squaring to further emphasize R-peaks
squared_data = diff_data .^ 2;

% Moving average integration
window_length = round(params.window_length * fs);
window = ones(1, window_length) / window_length;
integrated_data = conv(squared_data, window, 'same');

% Amplitude thresholding
threshold = params.threshold_ratio * max(integrated_data);

% Refractory period to avoid detecting multiple peaks within a short duration
refractory_half_period = round(params.refractory_period * fs / 2);


% Find the local peaks that satisfy both amplitude and width condition
peaks = zeros(1, length(data));
for k = 1 : length(peaks)
    search_win_start = max(1, k - refractory_half_period);
    search_win_stop = min(length(data), k + refractory_half_period);

    if integrated_data(k) > threshold && max(integrated_data(search_win_start : search_win_stop)) == integrated_data(k)
        peaks(k) = 1;
    end
end

% Get peak indexes
peak_indexes = find(peaks);

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
