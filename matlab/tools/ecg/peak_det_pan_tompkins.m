function [peaks, peak_indexes] = peak_det_pan_tompkins(ecg_data, fs, varargin)
% peak_det_pan_tompkins - R-peak detector based on Pan-Tompkins method.
% QRS Detection with Pan-Tompkins algorithm.
% Pan et al. 1985: A Real-Time QRS Detection Algorithm.
%   [peaks, peak_indexes] = peak_det_pan_tompkins(data, fs, varargin)
%
%   This function implements the Pan-Tompkins algorithm for R-peak detection
%   in ECG signals, with a simplified post-detection R-peak selection logic
%
%   Inputs:
%       data: Vector of input ECG data
%       fs: Sampling rate in Hz
%  (optional input):
%        - ecg_polarity: Optional. Search for positive (flag=1) or negative (flag=0) peaks.
%             By default, the skewness value of the bandpass-filtered signal determines the peak sign.
%        - qrs_width: expected maximal length of QRS complexes [ms]. Default 150ms.
%        - refracT:     refractory time for T-wave [ms] (default 360ms).
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Reference:
%       Pan J, Tompkins WJ. A real-time QRS detection algorithm. IEEE Trans
%       Biomed Eng. 1985;32(3):230-236. doi:10.1109/TBME.1985.325532
%
%   Reza Sameni, Sajjad Karimi 2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET
%
%
%   fc_low: BP filter lower cutoff frequency in Hz (default: 5.0 Hz)
%   fc_high: BP filter upper cutoff frequency in Hz (default: 15.0 Hz)

% Check optional input arguments
if nargin > 2 && ~isempty(varargin{1})
    ecg_polarity = varargin{1};
else
    ecg_polarity = [];
end

if nargin > 3 && ~isempty(varargin{2})
    qrs_width = varargin{2};
else
    qrs_width = 0.150;
end

if nargin > 4 && ~isempty(varargin{3})
    refracT = varargin{3};
else
    refracT = 0.360;
end



data = ecg_data(:);

fs_pt = 200;
total_delay = 0;

G = gcd( fs_pt,fs );
data = resample(data,fs_pt/G,fs/G);

%% Pan-Tomkins signal processing

% #1 remove the mean
data = data - mean(data);

% #2 low-pass filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
a = [1 -2 1];
filtered_data = filter(b,a,data);
total_delay = total_delay + 6;

% #3 high-pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
b = zeros(1,33);
b(1) = -1; b(17) = 32; b(33) = 1;
a = [1 1];
filtered_data = filter(b,a,filtered_data);    % Without Delay
total_delay = total_delay + 16;
filter_delay = total_delay;

% Differentiation to enhance R-peaks
% H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
b = [1 2 0 -2 -1].*(1/8)*fs;

diff_data = filter(b,1,filtered_data);
total_delay = total_delay + 2;

% Squaring to further emphasize R-peaks
squared_data = diff_data .^ 2;

% Moving average integration
window_length = round( qrs_width * fs_pt/2);
integrated_data = movmean(squared_data, [window_length,window_length]);

% remove delays & resample

integrated_data = [integrated_data(total_delay:end); repelem(integrated_data(end),total_delay-1,1)]';
integrated_data = resample(integrated_data,fs/G ,fs_pt/G);
integrated_data = integrated_data./std(integrated_data);

filtered_data = [filtered_data(filter_delay:end); repelem(filtered_data(end),filter_delay-1,1)]';
filtered_data = resample(filtered_data,fs/G ,fs_pt/G);

if isempty(ecg_polarity)
    filtered_data = filtered_data * sign(skew(filtered_data));
else
    filtered_data = filtered_data * sign(ecg_polarity-0.5);
end
filtered_data = filtered_data./std(filtered_data);


%% Adaptive thresholding

refracT = round( refracT * fs);

%=======================================================================
% ===================  learning phase 1 ==================================
%=======================================================================


% initializing thresholds based on first 3s of signal
temp_I = integrated_data(1:3*fs);
SPKI = max(temp_I);
NPKI = mean(temp_I);

thresholdI1 = NPKI + 0.25*(SPKI-NPKI);


temp_F = filtered_data(1:3*fs);
SPKF = max(temp_F);
NPKF = mean(temp_F);


T = length(integrated_data);
max_peak = T/fs * 3;

peaks = islocalmax(temp_I, 'MinSeparation',refracT, 'MinProminence',thresholdI1);

qrs_pos_init = find(peaks);

RR_init = median(diff(qrs_pos_init));

%=======================================================================
% =====================  detection phase ===============================
%=======================================================================

%initialization of detection phase
qrs_pos = nan(max_peak,1);
qrs_pos(1) = qrs_pos_init(1);

RR_recents = NaN(1,8);
RR_recents(1,1) = RR_init;


RR_selected = NaN(1,8);
RR_selected(1,1) = RR_init;


start_pos = qrs_pos(1,1) + find(temp_I(qrs_pos(1,1)+1:end)<thresholdI1,1,'first');
t = max(start_pos, qrs_pos(1,1)+refracT);

while t <= T-refracT/2

    RR1 = mean(RR_recents,'omitnan');
    RR2 = mean(RR_selected,'omitnan');

    %regulary heart rate check
    if (RR1 == RR2)
        factor = 1;
    else
        factor = 0.5;
    end


    %updating thresholds
    thresholdI1 = factor*(NPKI + (0.25*(SPKI-NPKI)));
    thresholdI2 = thresholdI1; %is automatically updated if necessary

    thresholdF1 = factor*(NPKF + (0.25*(SPKF-NPKF)));
    thresholdF2 = thresholdF1; %is automatically updated if necessary


    [t_I, above_thr, fall_thr] = segment_analysis(RR1, integrated_data, thresholdI2, start_pos, T );

    segment_normal = above_thr && fall_thr;
    segment_normal_I = segment_normal;

    %search back with threshold re-adjustment for possible R-peak
    thres_gain = 1;
    c = 0;
    while segment_normal==0 && c<5
        c=c+1;
        thres_gain = thres_gain*0.5;
        thresholdI2 = thresholdI2*thres_gain;
        [t_I, above_thr, fall_thr] = segment_analysis(RR1, integrated_data, thresholdI2, start_pos, T );
        segment_normal = above_thr && fall_thr;
    end

    [t_F, above_thr, fall_thr] = segment_analysis(RR1, filtered_data, thresholdF2, start_pos, T );
    segment_normal = above_thr && fall_thr;
    segment_normal_F = segment_normal;

    %search back with threshold re-adjustment for possible R-peak
    thres_gain = 1;
    c = 0;
    while segment_normal==0 && c<5
        c=c+1;
        thres_gain = thres_gain*0.5;
        thresholdF2 = thresholdF2*thres_gain;
        [t_F, above_thr, fall_thr] = segment_analysis(RR1, filtered_data, thresholdF2, start_pos, T );
        segment_normal = above_thr && fall_thr;
    end

    t = max(t_I, t_F);
    this_segment_index = start_pos:t;

    %begin new search from actual detection
    filt_data = filtered_data(this_segment_index);
    intg_data = integrated_data(this_segment_index);

    %exceeding both thresholds
    pos_common = intg_data>thresholdI2 & filt_data>thresholdF2;

    if  any(pos_common) % normal scenario

        pos_mask_I = find(pos_common);
        pos_mask_F = find(pos_common);

        neg_mask_I = find(~pos_common);
        neg_mask_F = find(~pos_common);

    else % for weak or irregular R-peak with back-search

        pos_I = intg_data>thresholdI2;
        pos_F = filt_data>thresholdF2;
        if ~any(pos_F)
            pos_F = pos_I;
        end

        pos_mask_I = find(pos_I);
        pos_mask_F = find(pos_F);

        neg_mask_I = find(~pos_I);
        neg_mask_F = find(~pos_F);

    end

    %find I-peaks
    idx_max_temp = islocalmax(intg_data(pos_mask_I) ,'MaxNumExtrema',1);
    if ~any(idx_max_temp)
        peakIs = max(intg_data(pos_mask_I));
    else
        peakIs = intg_data(pos_mask_I(idx_max_temp));
    end

    peakIn = mean(intg_data(neg_mask_I));


    %find R-peaks
    idx_max = islocalmax(filt_data(pos_mask_F) ,'MaxNumExtrema',1);
    idx_max = find(idx_max);
    if isempty(idx_max)
        [peakFs, idx_max] = max(filt_data(pos_mask_F));
    else
        peakFs = filt_data(pos_mask_F(idx_max));
    end

    peakFn = mean(intg_data(neg_mask_F));

    peakFs = max(0.01,peakFs);
    peakFn = max(0.01,peakFn);

    % update new R-peak position
    next_idx = find(isnan(qrs_pos),1,'first');
    qrs_pos(next_idx) = this_segment_index(pos_mask_F(idx_max));

    if next_idx==21
        yu=0;
    end

    if segment_normal_I && segment_normal_F
        SPKI = (0.125*peakIs)+(0.875*SPKI);
        NPKI = (0.125*peakIn)+(0.875*NPKI);
        SPKF = (0.125*peakFs)+(0.875*SPKF);
        NPKF = (0.125*peakFn)+(0.875*NPKF);
    else %double speed threshold adjustment after searchback
        SPKI = (0.25*peakIs)+(0.75*SPKI);
        NPKI = (0.25*peakIn)+(0.75*NPKI);
        SPKF = (0.25*peakFs)+(0.75*SPKF);
        NPKF = (0.25*peakFn)+(0.75*NPKF);
    end

    %check refractory/t-wave criterion and manipulate indices
    start_pos = t+1;
    start_pos = max(start_pos, qrs_pos(next_idx)+refracT);


    %calculate and sort in new RR-interval
    if next_idx >= 3
        RR_new = qrs_pos(next_idx)-qrs_pos(next_idx-1);

        RR_recents = circshift(RR_recents,1);
        RR_recents(1) = RR_new;

        if ((RR_new >= (0.92*RR1))&&(RR_new <= (1.16*RR1)))
            RR_selected = circshift(RR_selected,1);
            RR_selected(1) = RR_new;
        end

    end


    t = start_pos;

end

qrs_pos(isnan(qrs_pos))=[];
peaks = 0 * ecg_data;
peaks(qrs_pos) = 1;
peak_indexes = qrs_pos;

end


function [t, above_thr, fall_thr] = segment_analysis(RR_val, integrated_data, thresholdI2, start_pos, T )

counter_this = 0;
above_thr = 0;
fall_thr = 0;
t = start_pos;

while ((counter_this< 1.66 * RR_val) && ~(above_thr&&fall_thr)) && t < T
    counter_this = counter_this+1;
    t=t+1;

    if integrated_data(t)>thresholdI2
        above_thr = 1;
    end

    if integrated_data(t)<thresholdI2 &&  above_thr == 1
        fall_thr = 1;
    end

end

end
