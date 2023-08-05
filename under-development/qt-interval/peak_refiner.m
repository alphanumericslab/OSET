function [R_peak_process, R_quality_process, ecg_quality] = peak_refiner(ECG_raw, R_peak_indexes, fs)


% NOTCH FILTERING THE ECG
fc = 50.0; % powerline frequency
Qfactor = 45; % Q-factor of the notch filter
Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
[b,a] = iirnotch(Wo, BW); % design the notch filter
ECG_raw = filtfilt(b, a, ECG_raw); % zero-phase non-causal filtering

% ECG PRE-PROCESSING
fl_ECG = 10.0; % lower bandpass cut-off frequency for R-peaks
fh_ECG = 40.0; % upper bandpass cut-off frequency for R-peaks
ECG_hp = ECG_raw - lp_filter_zero_phase(ECG_raw, fl_ECG/fs);
ECG_bp = lp_filter_zero_phase(ECG_hp, fh_ECG/fs);
sigma = 5.0 * std(ECG_bp); % saturate peaks above k-sigma of the ECG
ECG_saturated = sigma * tanh(ECG_bp / sigma);



%% offline-processing
T = length(ECG_raw);
interval_width = round(fs/40);
L = length(R_peak_indexes);
R_peak_set = nan(L,2*interval_width+1);
R_peak_set_sat = nan(L,2*interval_width+1);

for r = 1:L

    this_peak_index = R_peak_indexes(r)-interval_width : R_peak_indexes(r)+interval_width;
    if this_peak_index(1)<1 || this_peak_index(end)>T
        continue;
    end
    R_peak_set(r,:) = ECG_raw(this_peak_index) - mean(ECG_raw(this_peak_index));
    R_peak_set_sat(r,:) = ECG_saturated(this_peak_index) - mean(ECG_saturated(this_peak_index));
end

pattern_R = median(R_peak_set,'omitnan');
pattern_R_sat = median(R_peak_set_sat,'omitnan');
R_quality = ((R_peak_set*pattern_R(:)/(2*interval_width+1)) ./ (std(pattern_R)*std(R_peak_set,[],2)) + (R_peak_set_sat*pattern_R_sat(:)/(2*interval_width+1)) ./ (std(pattern_R_sat)*std(R_peak_set_sat,[],2)))/2;

% first remove low-quality R-peaks
R_peak_high_quality = R_peak_indexes;

RR_intervals_ecg = diff(R_peak_high_quality); % RR-interval time series in samples
RR_intervals_ms = 1000*[RR_intervals_ecg(1);RR_intervals_ecg]/fs;

dif_for_missing = [0;diff(RR_intervals_ms)./movmedian(RR_intervals_ms(1:end-1),[10,10])];
dif_for_missing(dif_for_missing<0)=nan;

thr_val = 0.75;
index_missings = find(dif_for_missing>thr_val);

if ~isempty(index_missings)

    ecg_quality = "post-process";
    perc_thr = max(1,round(100*length(index_missings)/length(dif_for_missing)));
    thr = max(0.25, min(0.7,prctile(R_quality,perc_thr) ));
    R_peak_high_quality(R_quality<thr)=[];

    % fine-tune high-quality peaks with ECG_raw
    for r = 1:length(R_peak_high_quality)
        this_peak_index = R_peak_high_quality(r)-interval_width : R_peak_high_quality(r)+interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ECG_raw(this_peak_index) - mean(ECG_raw(this_peak_index));

        [cxy,lags] = xcov(temp_ecg,pattern_R,interval_width,'biased');
        [~,idx_max]= max(cxy);
        %             [~,idx_max]= max(temp_ecg*sign(pattern_R(interval_width+1)));
        %             lags = -interval_width : interval_width;
        R_peak_high_quality(r) =  R_peak_high_quality(r) + lags(idx_max);
    end

    RR_intervals_ecg = diff(R_peak_high_quality); % RR-interval time series in samples
    RR_intervals_ms = 1000*[RR_intervals_ecg(1);RR_intervals_ecg]/fs;
    med_RR_intervals_ms = movmedian(RR_intervals_ms,[10,10]);
    dif_for_missing = [0;diff(RR_intervals_ms)./med_RR_intervals_ms(1:end-1)];
    dif_for_missing(dif_for_missing<0)=nan;

    thr_val = 0.75;
    index_missings = find(dif_for_missing>thr_val);
    num_missing = round(dif_for_missing(index_missings));

    R_peak_process = R_peak_high_quality;
    R_peak_insert = [];
    for r = 1:length(index_missings)
        imput_index = index_missings(r);
        expected_sample = round((fs/1000)*med_RR_intervals_ms(imput_index));
        R_peak_insert = [ R_peak_insert; R_peak_process(imput_index-1)+(1:num_missing(r))'*expected_sample] ;
    end


    % fine tune this peaks with ECG_raw
    for r = 1:length(R_peak_insert)
        this_peak_index = R_peak_insert(r)-1*interval_width : R_peak_insert(r)+1*interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ECG_raw(this_peak_index) - mean(ECG_raw(this_peak_index));

        [cxy,lags] = xcov(temp_ecg,pattern_R,interval_width,'biased');
        [~,idx_max]= max(cxy);
        %             [~,idx_max]= max(temp_ecg*sign(pattern_R(interval_width+1)));
        %             lags = -1*interval_width : 1*interval_width;
        R_peak_insert(r) =  R_peak_insert(r) + lags(idx_max);
    end

    R_peak_process = sort(unique([R_peak_process;R_peak_insert]),'ascend');

    %% second time
    R_peak_high_quality = R_peak_process;
    RR_intervals_ecg = diff(R_peak_high_quality); % RR-interval time series in samples
    RR_intervals_ms = 1000*[RR_intervals_ecg(1);RR_intervals_ecg]/fs;
    med_RR_intervals_ms = movmedian(RR_intervals_ms,[10,10]);

    dif_for_missing = [0;diff(RR_intervals_ms)./med_RR_intervals_ms(1:end-1)];
    dif_for_missing(dif_for_missing<0)=nan;

    thr_val =0.75;
    index_missings = find(dif_for_missing>thr_val);
    num_missing = round(dif_for_missing(index_missings));

    R_peak_process = R_peak_high_quality;
    R_peak_insert = [];
    for r = 1:length(index_missings)
        imput_index = index_missings(r);
        expected_sample = round((fs/1000)*med_RR_intervals_ms(imput_index));
        R_peak_insert = [ R_peak_insert; R_peak_process(imput_index-1)+(1:num_missing(r))'*expected_sample] ;
    end


    % fine tune this peaks with ECG_raw
    for r = 1:length(R_peak_insert)
        this_peak_index = R_peak_insert(r)-1*interval_width : R_peak_insert(r)+1*interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ECG_raw(this_peak_index) - mean(ECG_raw(this_peak_index));

        [cxy,lags] = xcov(temp_ecg,pattern_R,interval_width,'biased');
        [~,idx_max]= max(cxy);
        %             [~,idx_max]= max(temp_ecg*sign(pattern_R(interval_width+1)));
        %             lags = -1*interval_width : 1*interval_width;
        R_peak_insert(r) =  R_peak_insert(r) + lags(idx_max);
    end

    R_peak_process = sort(unique([R_peak_process;R_peak_insert]),'ascend');

    L = length(R_peak_process);
    R_peak_set = nan(L,2*interval_width+1);

    for r = 1:L

        this_peak_index = R_peak_process(r)-interval_width : R_peak_process(r)+interval_width;
        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end
        R_peak_set(r,:) = ECG_raw(this_peak_index) - mean(ECG_raw(this_peak_index));

    end


    pattern_R = median(R_peak_set,'omitnan');
    R_quality_process = (R_peak_set*pattern_R(:)/(2*interval_width+1)) ./ (std(pattern_R)*std(R_peak_set,[],2));

else
    ecg_quality = "clean";
    % fine-tune high-quality peaks with ECG_raw
    for r = 1:length(R_peak_high_quality)
        this_peak_index = R_peak_high_quality(r)-interval_width : R_peak_high_quality(r)+interval_width;

        if this_peak_index(1)<1 || this_peak_index(end)>T
            continue;
        end

        temp_ecg = ECG_raw(this_peak_index) - mean(ECG_raw(this_peak_index));

        [~,idx_max]= max(temp_ecg*sign(pattern_R(interval_width+1)));
        lags = -1*interval_width : 1*interval_width;
        R_peak_high_quality(r) =  R_peak_high_quality(r) + lags(idx_max);
    end
    R_peak_process = R_peak_high_quality;
    R_quality_process = R_quality;

end

