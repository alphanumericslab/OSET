function rr_intervals_ecg = em_interval_calc(ecg_rpeaks_index, thr_max, thr_fast_max,max_lag)

if nargin<2
    thr_max = 0.5;
    thr_fast_max = 0.4;
    max_lag= 20;
elseif nargin<3
    thr_fast_max = 0.4;
    max_lag= 20;
end

if size(ecg_rpeaks_index,2)==1
    rr_intervals_ecg = diff(ecg_rpeaks_index); % RR-interval time series in samples
elseif size(ecg_rpeaks_index,2)==2
    rr_intervals_ecg = ecg_rpeaks_index(:,2)-ecg_rpeaks_index(:,1);
end

md_rr_intervals_ecg = movmedian(rr_intervals_ecg,[max_lag,max_lag],'omitnan');
diff_outliers = rr_intervals_ecg./md_rr_intervals_ecg-1;
ind_nan_rr = find(abs(diff_outliers)>thr_max);
rr_intervals_ecg(ind_nan_rr) = md_rr_intervals_ecg(ind_nan_rr);

if ~isempty(ind_nan_rr)
    md_rr_intervals_ecg = movmedian(rr_intervals_ecg,[ceil(max_lag/3),ceil(max_lag/3)],'omitnan');
    diff_outliers = rr_intervals_ecg./md_rr_intervals_ecg-1;
    ind_nan_rr = find(abs(diff_outliers)>thr_fast_max);
    rr_intervals_ecg(ind_nan_rr) = md_rr_intervals_ecg(ind_nan_rr);

    diff_outliers = rr_intervals_ecg-md_rr_intervals_ecg;
    std_rr_intervals_ecg = movstd(diff_outliers,[100,100],'omitnan');
    z_score_rr = abs(diff_outliers./std_rr_intervals_ecg);
    ind_nan_zrr = find(abs(z_score_rr)>5);
    rr_intervals_ecg(ind_nan_zrr) = md_rr_intervals_ecg(ind_nan_zrr);
end

 rr_intervals_ecg = fillmissing(rr_intervals_ecg,'linear');
