function [peaks, peak_indexes, width] = peak_det_modified_pan_tompkins(data, fs, varargin)
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
%       ksigma: Optional. Saturates peaks at ksigma x STD of the signal after energy envelope calculation (default = 12)
%       flag: Optional. Search for positive (flag=1) or negative (flag=0) peaks.
%             By default, the maximum absolute value of the signal determines the peak sign.
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%       width: Rise to fall width of the signal's peak bump (in samples)
%
%   Revision History:
%       2006: First release
%       2023: Added ksigma saturation feature, removed a bug in baseline
%           wander remover filter and renamed from deprecated version
%           PeakDetection2
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
    wlen = 0.150; % moving average window length 150ms
end

if nargin > 3 && ~isempty(varargin{2})
    fp1 = varargin{2};
else
    fp1 = 10; % First zero of the HP filter is placed at f = 10Hz;
end

if nargin > 4 && ~isempty(varargin{3})
    fp2 = varargin{3};
else
    fp2 = 33.3; % First zero of the LP filter is placed at f = 33.3Hz;
end

if nargin > 5 && ~isempty(varargin{4})
    th = varargin{4};
else
    th = 0.2;
end

if nargin > 6 && ~isempty(varargin{5})
    ksigma = varargin{5};
else
    ksigma = 12;
end

if nargin > 7 && ~isempty(varargin{6})
    flag = varargin{6};
else
    flag = abs(max(data)) > abs(min(data));
end

N = length(data);
data = data(:)';

L1 = round(fs / fp2); 
L2 = round(fs / fp1); 

x0 = data - lp_filter_zero_phase(data, 0.05 / fs);

% LP filter
x = filter([1, zeros(1, L1 - 1), -1], [L1, -L1], x0);
x = filter([1, zeros(1, L1 - 1), -1], [L1, -L1], x);
x = [x(L1:end), zeros(1, L1 - 1) + x(end)]; % Lag compensation

% HP filter
y = filter([L2 - 1, -L2, zeros(1, L2 - 2), 1], [L2, -L2], x);

% Differentiation
z = diff([y(1) y]);

% Squaring
w = z .^ 2;

% Moving average
L3 = round(fs * wlen);
v = filter([1, zeros(1, L3 - 1), -1], [L3, -L3], w);
v = [v(round(L3 / 2):end), zeros(1, round(L3 / 2) - 1) + v(end)]; % Group-delay lag compensation

v_sat = tanh_saturation(v, ksigma); % Saturate peaks

p = v_sat > (th * max(v_sat)); % Potential peaks

% Edge detection
rising = find(diff([0, p]) == 1);      % Rising edges
falling = find(diff([p, 0]) == -1);    % Falling edges

if length(rising) == length(falling) - 1
    rising = [1; rising];
elseif length(rising) == length(falling) + 1
    falling = [falling; N];
end

peak_indexes = zeros(length(rising), 1);
width = zeros(length(rising), 1);

if flag
    for i = 1:length(rising)
        [~, mx] = max(data(rising(i) : falling(i)));
        peak_indexes(i) = mx - 1 + rising(i);
        width(i) = falling(i) - rising(i);
    end
else
    for i = 1:length(rising)
        [~, mn] = min(data(rising(i) : falling(i)));
        peak_indexes(i) = mn - 1 + rising(i);
        width(i) = falling(i) - rising(i);
    end
end

peaks = zeros(1, N);
peaks(peak_indexes) = 1;
end
