function [peaks, peak_indexes, R_prominence] = peak_det_advanced_pan_tompkins(data, fs, varargin)
%
% peak_det_modified_pan_tompkins - R-peak detector based on modified
%   Pan-Tompkins method. The filters and post-detection R-peak selection
%   logic differ from the original algorithm
%
%   [peaks, peak_indexes, width] = peak_det_modified_pan_tompkins(data, fs, wlen, fp1, fp2, th, ksigma, flag)
%
%   Inputs:
%       data: Vector of input data
%       fs: Sampling rate
%       wlen: Optional. Moving average window length (default = 150ms)
%       fp1: Optional. Lower cut-off frequency (default = 10Hz)
%       fp2: Optional. Upper cut-off frequency (default = 33.3Hz)
%       th: Optional. Detection threshold (default = 0.2)
%       ksigma: Optional. Saturates peaks at ksigma x STD of the signal after energy envelope calculation (default = 5)
%       flag: Optional. Search for positive (flag=1) or negative (flag=0) peaks.
%             By default, the maximum absolute value of the signal determines the peak sign.
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Revision History:
%       2023: First release
%
%   Reference:
%       Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
%       Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532
%
%   Reza Sameni, 2008-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Set default values for optional arguments
if nargin > 2 && ~isempty(varargin{1})
    wlen = varargin{1};
else
    wlen = 0.02; % moving average window length 150ms
end

if nargin > 3 && ~isempty(varargin{2})
    fp1 = varargin{2};
else
    fp1 = 10; % strat-band for band-pass filter f = 10Hz;
end

if nargin > 4 && ~isempty(varargin{3})
    fp2 = varargin{3};
else
    fp2 = 33.3; % stop-band for band-pass filter f = 33.3Hz;
end

if nargin > 5 && ~isempty(varargin{4})
    th = varargin{4};
else
    th = 0.1;
end

if nargin > 6 && ~isempty(varargin{6})
    ksigma = varargin{6};
else
    ksigma = 5;
end

if nargin > 7 && ~isempty(varargin{5})
    flag = varargin{5};
else
    flag = -1;
end



% N = length(data);
data = data(:)';

ecg_base = data - lp_filter_zero_phase(data, 2/fs);
ecg_base = tanh_saturation(ecg_base, 3); % Saturate v_max

x = data - lp_filter_zero_phase(data, fp1/fs);
y = lp_filter_zero_phase(x, fp2/fs);
sigma = ksigma * std(y); % saturate peaks above k-sigma of the ECG
y = sigma * tanh(y / sigma);

% Differentiation
d_order = 1;
z = abs(diff([repmat(y(1),1,d_order) y],d_order));


% Squaring
w = z .^ 2;

% Moving average
L3 = ceil(fs / 50);

v = movmax(w,[L3 ,L3]);

% Another Differentiation
v = abs(diff([repmat(v(1),1,d_order) v],d_order));

v = movsum(v,[L3+1 ,L3+1]);

v = movmax(v,[L3+1 ,L3+1]);
v = movmax(v,[L3+1 ,L3+1]);


v_sat = tanh_saturation(v, ksigma); % Saturate v_max
v_sat = tanh_saturation(v_sat, ksigma); % Saturate v_max
v_sat = tanh_saturation(v_sat, ksigma); % Saturate v_max

data_skew = v_sat .* y/max(abs(y));

if flag==-1
    if(skew(data_skew) > 0) % find the ECG polarity based on skewness sign
        ecg_polarity = 1;
    else
        ecg_polarity = -1;
    end
else
    if flag ==1
        ecg_polarity = 1;
    else
        ecg_polarity = -1;
    end

end

data_rpeak_enhanced = v_sat .* ecg_base/max(abs(ecg_base));

[peaks, R_prominence] = islocalmax(ecg_polarity*data_rpeak_enhanced, 'MinSeparation',15*L3, 'MinProminence',(th * max(v_sat)));
peak_indexes = find(peaks);
R_prominence = R_prominence(peak_indexes)/max(v_sat);


end
