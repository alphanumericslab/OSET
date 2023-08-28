function [peak_refined_indexes, peak_quality_scores, ecg_quality] = peak_refiner_matched_filter(ecg_raw, peak_indexes, fs, params)

% peak_refiner - R-peak refiner for an existing R-peak detector method with R-peaks in R_peak_indexes
% Syntax:
%   [R_peak_refined, R_peak_quality, ecg_quality] = peak_refiner(ECG_raw, R_peak_indexes, fs)
%
% Method:
%   refine R-peaks in the ECG signal for an existing R-peak detector method
%   using matched filter concept & RR series analysis. 
%   The function returns corresponding refined indexes (peak_refined_indexes).
%
% Inputs:
%   ecg_raw: Vector of the input ECG signal.
%   peak_indexes: Indexes of the detected R-peaks in the input signal.
%   fs: Sampling frequency in Hz.
%   params.f_powerline : powerline frequency in Hz (default = 60Hz).

% Outputs:
%   peak_refined_indexes: Indexes of the detected R-peaks in the input signal.
%   peak_quality_scores: Score between [-1,1] that indicates the quality of detected peaks
%   ecg_quality: Indicate the quality of extracted peak_indexes for the input ECG

% Revision History:
% 2023: First release.
%
% Reza Sameni, Sajjad Karimi 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET
if nargin>3
    if isfield(params,'f_powerline')
        fc = params.f_powerline;
    else
        fc = 60.0; % powerline frequency
    end

else
    fc = 60.0; % powerline frequency
end

% Notch filtering the ECG
Qfactor = 45; % Q-factor of the notch filter
Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
[b,a] = iirnotch(Wo, BW); % design the notch filter
ecg_raw = filtfilt(b, a, ecg_raw); % zero-phase non-causal filtering

% ECG PRE-PROCESSING
fl_ECG = 10.0; % lower bandpass cut-off frequency for R-peaks
fh_ECG = 40.0; % upper bandpass cut-off frequency for R-peaks
ECG_hp = ecg_raw - lp_filter_zero_phase(ecg_raw, fl_ECG/fs);
ECG_bp = lp_filter_zero_phase(ECG_hp, fh_ECG/fs);
sigma = 5.0 * std(ECG_bp); % saturate peaks above k-sigma of the ECG
ECG_saturated = sigma * tanh(ECG_bp / sigma);


%% Start processing R-peaks

T = length(ecg_raw);
interval_width = round(fs/40);
L = length(peak_indexes);
R_peak_set = nan(L,2*interval_width+1);
R_peak_set_sat = nan(L,2*interval_width+1);

for r = 1:L
    this_peak_index = peak_indexes(r)-interval_width : peak_indexes(r)+interval_width;
    if this_peak_index(1)<1 || this_peak_index(end)>T
        continue;
    end
    R_peak_set(r,:) = ecg_raw(this_peak_index) - mean(ecg_raw(this_peak_index));
    R_peak_set_sat(r,:) = ECG_saturated(this_peak_index) - mean(ECG_saturated(this_peak_index));
end

pattern_R = median(R_peak_set,'omitnan');
pattern_R_sat = median(R_peak_set_sat,'omitnan');
peak_quality_scores = ((R_peak_set*pattern_R(:)/(2*interval_width+1)) ./ (std(pattern_R)*std(R_peak_set,[],2)) + (R_peak_set_sat*pattern_R_sat(:)/(2*interval_width+1)) ./ (std(pattern_R_sat)*std(R_peak_set_sat,[],2)))/2;

% first remove low-quality R-peaks
R_peak_high_quality = peak_indexes;

RR_intervals_ecg = diff(R_peak_high_quality); % RR-interval time series in samples
RR_intervals_ms = 1000*[RR_intervals_ecg(1);RR_intervals_ecg]/fs;

diff_for_missing = [0;diff(RR_intervals_ms)./movmedian(RR_intervals_ms(1:end-1),[10,10])];
diff_for_missing(diff_for_missing<0)=nan;

thr_val = 0.75;
index_missings = find(diff_for_missing>thr_val);

% Check ecg_quality
% Decide to start post-processing or just fine-tune peaks with a max-operator for clean data

if ~isempty(index_missings)

    % Start post-processing of R-peaks
    ecg_quality = "post-process";
    pattern_R_extend = [nan(1,interval_width), pattern_R, nan(1,interval_width)];

    perc_thr = max(1,round(100*length(index_missings)/length(diff_for_missing)));
    thr = max(0.25, min(0.7,prctile(peak_quality_scores,perc_thr) ));
    R_peak_high_quality(peak_quality_scores<thr)=[];

    % fine-tune high-quality peaks with matched filter
    for r = 1:length(R_peak_high_quality)
        this_peak_index = R_peak_high_quality(r)-interval_width : R_peak_high_quality(r)+interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end
        temp_ecg = ecg_raw(this_peak_index) - mean(ecg_raw(this_peak_index));

        % Matched filter
        [cxy, lags] = xcov(temp_ecg, pattern_R, interval_width,'biased');
        [~,idx_max]= max(cxy);
        R_peak_high_quality(r) =  R_peak_high_quality(r) + lags(idx_max);
    end

    RR_intervals_ecg = diff(R_peak_high_quality); % RR-interval time series in samples
    RR_intervals_ms = 1000*[RR_intervals_ecg(1);RR_intervals_ecg]/fs;
    med_RR_intervals_ms = movmedian(RR_intervals_ms,[10,10]);
    diff_for_missing = [0;diff(RR_intervals_ms)./med_RR_intervals_ms(1:end-1)];
    diff_for_missing(diff_for_missing<0)=nan;

    thr_val = 0.75;
    index_missings = find(diff_for_missing>thr_val);
    num_missing = round(diff_for_missing(index_missings));

    peak_refined_indexes = R_peak_high_quality;
    R_peak_inserted = [];
    for r = 1:length(index_missings)
        imput_index = index_missings(r);
        expected_sample = round((fs/1000)*med_RR_intervals_ms(imput_index));
        R_peak_inserted = cat(1, R_peak_inserted, peak_refined_indexes(imput_index-1)+(1:num_missing(r))'*expected_sample) ;
    end


    % fine-tune new inserted peaks using matched filter and raw ECG
    for r = 1:length(R_peak_inserted)
        this_peak_index = R_peak_inserted(r)-2*interval_width : R_peak_inserted(r)+2*interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ecg_raw(this_peak_index) - mean(ecg_raw(this_peak_index));

        % Matched filter
        %         [cxy, lags] = xcov(temp_ecg, pattern_R_extend, interval_width,'unbiased');
        cxy = autocrosscorrelation_missing( temp_ecg, pattern_R_extend, round(1.5*interval_width));
        lags = -round(1.5*interval_width) : round(1.5*interval_width);
        [val_max,idx_max]= max(cxy);
        if val_max>0.3
            R_peak_inserted(r) =  R_peak_inserted(r) + lags(idx_max);
        end
    end

    peak_refined_indexes = sort(unique([peak_refined_indexes;R_peak_inserted]),'ascend');

    % Start post-processing of R-peaks one more time
    R_peak_high_quality = peak_refined_indexes;
    RR_intervals_ecg = diff(R_peak_high_quality); % RR-interval time series in samples
    RR_intervals_ms = 1000*[RR_intervals_ecg(1);RR_intervals_ecg]/fs;
    med_RR_intervals_ms = movmedian(RR_intervals_ms,[10,10]);

    diff_for_missing = [0;diff(RR_intervals_ms)./med_RR_intervals_ms(1:end-1)];
    diff_for_missing(diff_for_missing<0)=nan;

    thr_val =0.75;
    index_missings = find(diff_for_missing>thr_val);
    num_missing = round(diff_for_missing(index_missings));

    peak_refined_indexes = R_peak_high_quality;
    R_peak_inserted = [];
    for r = 1:length(index_missings)
        imput_index = index_missings(r);
        expected_sample = round((fs/1000)*med_RR_intervals_ms(imput_index));
        R_peak_inserted = cat(1, R_peak_inserted, peak_refined_indexes(imput_index-1)+(1:num_missing(r))'*expected_sample) ;
    end


    % fine-tune this peaks with ECG_raw
    for r = 1:length(R_peak_inserted)
        this_peak_index = R_peak_inserted(r)-2*interval_width : R_peak_inserted(r)+2*interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ecg_raw(this_peak_index) - mean(ecg_raw(this_peak_index));

        % Matched filter
        %         [cxy, lags] = xcov(temp_ecg, pattern_R_extend, interval_width);
        cxy = autocrosscorrelation_missing( temp_ecg, pattern_R_extend, round(1.5*interval_width));
        lags = -round(1.5*interval_width) : round(1.5*interval_width);
        [val_max,idx_max]= max(cxy);
        if val_max>0.3
            R_peak_inserted(r) =  R_peak_inserted(r) + lags(idx_max);
        end

    end

    peak_refined_indexes = sort(unique([peak_refined_indexes;R_peak_inserted]),'ascend');

    L = length(peak_refined_indexes);
    R_peak_set = nan(L,2*interval_width+1);

    for r = 1:L

        this_peak_index = peak_refined_indexes(r)-interval_width : peak_refined_indexes(r)+interval_width;
        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end
        R_peak_set(r,:) = ecg_raw(this_peak_index) - mean(ecg_raw(this_peak_index));
    end

    pattern_R = median(R_peak_set,'omitnan');
    peak_quality_scores = (R_peak_set*pattern_R(:)/(2*interval_width+1)) ./ (std(pattern_R)*std(R_peak_set,[],2));

else

    ecg_quality = "clean";
    % fine-tune high-quality peaks with ECG_raw
    for r = 1:length(R_peak_high_quality)
        this_peak_index = R_peak_high_quality(r)-interval_width : R_peak_high_quality(r)+interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ecg_raw(this_peak_index) - mean(ecg_raw(this_peak_index));

        [~,idx_max]= max(temp_ecg*sign(pattern_R(interval_width+1)));
        lags = -1*interval_width : 1*interval_width;
        R_peak_high_quality(r) =  R_peak_high_quality(r) + lags(idx_max);
    end
    peak_refined_indexes = R_peak_high_quality;
    %     peak_quality_scores = peak_quality_scores;

end



function Cxy = autocrosscorrelation_missing( X , Y , Lags)

%Cxy(T)=E{X(t)Y(t-T)}


Cxy = zeros(1,2*Lags+1);
Cxy = Cxy *nan;
Cxy(Lags+1)=1;

for k = 1:Lags

    temp1 = X(1:end-k);
    temp2= Y(k+1:end);

    ind_nan1 = isnan(temp1);
    ind_nan2 = isnan(temp2);
    indnan = ind_nan1|ind_nan2;

    temp1 = temp1(~indnan);
    temp2 = temp2(~indnan);

    if(length(temp1)>2)
        Cxy(Lags+1-k) = (temp1(:)'-mean(temp1(:)))*(temp2(:)-mean(temp2(:)))/sqrt((length(temp1)-1)^2*var(temp1)*var(temp2));
        %Cxy(Lags+1-k) = (temp1(:)'-mean(temp1(:)))*(temp2(:)-mean(temp2(:)));
    end

end


for k = 0:Lags

    temp1 = X(k+1:end);
    temp2= Y(1:end-k);

    ind_nan1 = isnan(temp1);
    ind_nan2 = isnan(temp2);
    indnan = ind_nan1|ind_nan2;

    temp1 = temp1(~indnan);
    temp2 = temp2(~indnan);

    if(length(temp1)>2)
        Cxy(k+1+Lags) =  (temp1(:)'-mean(temp1(:)))*(temp2(:)-mean(temp2(:)))/sqrt(length(temp1)^2*var(temp1)*var(temp2));
        %Cxy(k+1+Lags) =  (temp1(:)'-mean(temp1(:)))*(temp2(:)-mean(temp2(:)));
    end

end

