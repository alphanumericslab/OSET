function [score1, score2, peaks1, peaks2, rank1, rank2] = sqi_multi_matched_filter(x,f,fs, varargin)
% sqi_multi_matched_filter - Signal quality index based on three fixed
%   template matched filter outputs scatter and average waveform standard
%   deviation
%
% Usage:
%   [score1, score2, peaks1, peaks2, rank1, rank2] = sqi_multi_matched_filter(x, f, fs, num_peak_detection_itr, template_type, template_params, plot_templates)
%
% Inputs:
%   x: Input data array (channels x samples).
%   f: Frequency vector.
%   fs: Sampling frequency.
%   num_peak_detection_itr (optional): Number of iterations for peak
%       detection (default: 1). 
%   template_type (optional): Template type (default: 0).
%       template_type 0: parameters preselected for fetal ECG
%       template_type 1: parameters preselected for adult ECG
%       template_type 2: provide QRS effective width and STD in template_params
%   template_params (optional): bi-vector in the form
%       [QRS_wave_effective_width, QRS_std], where QRS_wave_effective_width is
%       the total length of the QRS template in seconds and QRS_std is the
%       standard deviation of the wave
%   plot_templates (optional): Plot templates flag (default: 0).
%
% Outputs:
%   score1: Quality score based on amplitude ratio of spikes over non-spikes.
%   score2: Quality score based on variance of small beat variations.
%   peaks1: R-peaks from the best quality channel according to score1.
%   peaks2: R-peaks from the best quality channel according to score2.
%   rank1: Channel ranking based on score1.
%   rank2: Channel ranking based on score2.
%
% Method: The sqi_multi_matched_filter method computes a signal quality
%   index by applying template matched filtering to multichannel data. It
%   uses predefined templates to filter the data and extract relevant
%   features related to the R-peaks of the signals. The filtered outputs are
%   combined to calculate a power envelope, and local peaks are detected.
%   Quality scores are then computed based on the amplitude ratio of spikes
%   and non-spikes and the variance of small beat variations around the
%   R-peaks. The method ranks the channels based on these scores and extracts
%   R-peaks from the best quality channels. See reference.
% 
% References:
%   Biglari, H., & Sameni, R. (2016). Fetal motion estimation from
%   noninvasive cardiac signal recordings. In Physiological Measurement (Vol.
%   37, Issue 11, pp. 2003–2023). IOP Publishing.
%   https://doi.org/10.1088/0967-3334/37/11/2003
%
% Revision History:
%   2015: First release
%   2023: Renamed from deprecated version ChannelIndex12
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Handle optional arguments
if nargin > 3 && ~isempty(varargin{1})
    num_peak_detection_itr = varargin{1};
else
    num_peak_detection_itr = 1;
end

if nargin > 4 && ~isempty(varargin{2})
    template_type = varargin{2};
else
    template_type = 0;
end

if nargin > 5 && template_type~= 0 && ~isempty(varargin{3})
    param = varargin{3};
    beat_onset = -param(1)/2; % s
    beat_offset = param(1)/2; % s
    beat_std = param(2); % s
else
    beat_onset = -0.03; % s
    beat_offset = 0.03; % s
    beat_std = 0.01; % s
end

if nargin > 6 && ~isempty(varargin{4})
    plot_templates = varargin{4};
else
    plot_templates = 0;
end

% Define templates based on template_type
switch template_type
    case 0
        t1 = -.02 : 1/fs : .02; %s
        template1 = exp(-t1.^2/(0.007)^2);
        
        t2 = t1(2 : end) - 1/(2*fs); %s
        template2 = diff(exp(-t1.^2/(0.007)^2));
        
        t3 = t1; %s
        template3 = -.2*exp(-(t1+.007).^2/(0.005)^2) + exp(-t1.^2/(0.005)^2) + -.2*exp(-(t1-.007).^2/(0.005)^2);
    case 1
        t1 = beat_onset : 1/fs : beat_offset;
        template1 = exp(-t1 .^ 2 / beat_std ^ 2);
        
        t2 = t1(2:end) - 1/(2*fs);
        template2 = diff(template1);
        
        t3 = t1(3:end) - 1/(fs);
        template3 = diff(template1, 2);
end

% Plot templates if plot_templates flag is set
if plot_templates
    figure;
    plot(t1,template1);
    hold on
    plot(t2,template2);
    plot(t3,template3);
    grid;
    legend('template #1', 'template #2', 'template #3')
    title('R-peak detection matched filter templates')
end

L = size(x,1);
N = size(x,2);

score1 = zeros(1,L);
score2 = zeros(1,L);
ppeaks = zeros(L, N);
for i = 1:L
    xx = x(i,:); % select channel

    % First template matched filter
    L1 = length(template1);
    h1 = template1(end:-1:1);
    w1 = floor(L1/2);
    r1 = filter(h1,sqrt(sum(h1.^2)*sum(xx.^2)),[xx zeros(1,w1-1)]);
    r1 = r1(w1:N+w1-1);

    % Second template matched filter
    L2 = length(template2);
    h2 = template2(end:-1:1);
    w2 = floor(L2/2);
    r2 = filter(h2,sqrt(sum(h2.^2)*sum(xx.^2)),[xx zeros(1,w2-1)]);
    r2 = r2(w2:N+w2-1);

    % Third template matched filter
    L3 = length(template3);
    h3 = template3(end:-1:1);
    w3 = floor(L3/2);
    r3 = filter(h3,sqrt(sum(h3.^2)*sum(xx.^2)),[xx zeros(1,w3-1)]);
    r3 = r3(w3:N+w3-1);

    % 3D matched filter output power envelope
    r = sqrt(r1.^2 + r2.^2 + r3.^2);

    % Find the local peaks of the power envelope
    ppeaks(i, :) = peak_detection_local_search(r, f/fs, 1, num_peak_detection_itr);

    % Combine all nearby peaks to form a pulse corresponding to the R-peak of length LL (widest template width)
    LL = max([L1 L2 L3]);
    peakpulses = conv(ppeaks(i, :), ones(1, LL), 'same');
    I_spikes = peakpulses == 1;
    I_nonspikes = peakpulses == 0;

    % First score: amplitude ratio of spikes over non-spikes
    score1(i) = median(r(I_spikes))/median(r(I_nonspikes));

    % Second score: small beat variances around the R-peaks
    I = find(ppeaks(i, :) == 1);
    I(1) = []; I(end) = []; % Exclude first and last beats to avoid boundary effects
    blocks = zeros(length(I), 2*LL);
    for k = 1 : length(I)
        blocks(k, :) = r(I(k) - LL : I(k) + LL - 1)/std(r);
    end
    sample_std = std(blocks, [], 1);
    score2(i) = 1/mean(sample_std);
end

% Rank channels based on scores
[~, rank1] = sort(score1, 'descend');
peaks1 = ppeaks(rank1(1), :);

[~, rank2] = sort(score2, 'descend');
peaks2 = ppeaks(rank2(1), :);
