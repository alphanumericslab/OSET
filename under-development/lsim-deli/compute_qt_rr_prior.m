function qt_rr_prior = compute_qt_rr_prior(data_bias, ecg_rpeaks_index, ...
    ecg_qrson_index, ecg_qrsoff_index, ...
    beat_quality_twave, avg_intervals_ecg, rr_intervals_ecg, ...
    index_clustering, num_cls, fs, win_sample_T, quality_thr, rr_percentile)
%COMPUTE_QT_RR_PRIOR  Estimate QT and RR from clean high-quality beats per cluster
%
%   qt_rr_prior = compute_qt_rr_prior(data_bias, ecg_rpeaks_index,
%       ecg_qrson_index, ecg_qrsoff_index, beat_quality_twave,
%       avg_intervals_ecg, rr_intervals_ecg, index_clustering,
%       num_cls, fs, win_sample_T, quality_thr, rr_percentile)
%
%   Selects a high-quality subset of beats from each morphology cluster,
%   preferring beats with lower heart rate (longer RR intervals), and
%   estimates QT and RR intervals. The output serves as a prior to adjust
%   the qt_curve for robust T-wave offset detection during noisy or
%   exercise periods.
%
%   Inputs:
%       data_bias          - Bias-corrected ECG signal (1 x T)
%       ecg_rpeaks_index   - R-peak sample positions (N x 1)
%       ecg_qrson_index    - QRS onset positions (N x 1), i.e. positions.QRSon
%       ecg_qrsoff_index   - QRS offset positions (N x 1), i.e. positions.QRSoff
%       beat_quality_twave - Beat quality score per beat (N x 1), higher is better
%       avg_intervals_ecg  - Smoothed RR intervals in samples (N x 1)
%       rr_intervals_ecg   - RR intervals in samples (N x 1)
%       index_clustering   - (M x 2) matrix [inner_beat_index, cluster_id]
%                            inner_beat_index is 1-based for beats 2..N-1
%       num_cls            - Unique cluster IDs vector
%       fs                 - Sampling rate (Hz)
%       win_sample_T       - Window size for T-wave features (samples)
%       quality_thr        - (optional) Min quality for clean subset, default 0.75
%       rr_percentile      - (optional) Percentile of avg RR above which beats
%                            are selected (higher = prefer lower HR), default 50
%
%   Output:
%       qt_rr_prior - Structure with per-cluster results:
%         .cluster(c).QT              - QT intervals (samples) for clean beats
%         .cluster(c).RR              - RR intervals (samples) for clean beats
%         .cluster(c).median_QT       - Median QT (samples)
%         .cluster(c).median_RR       - Median RR (samples)
%         .cluster(c).qt_template     - Template-based QT estimate (samples)
%         .cluster(c).num_beats       - Number of clean beats used
%         .cluster(c).cluster_id      - Cluster identifier
%         .cluster(c).clean_beat_idx  - Beat indices in ecg_rpeaks_index
%
%   Sajjad Karimi 2024

if nargin < 12 || isempty(quality_thr)
    quality_thr = 0.75;
end
if nargin < 13 || isempty(rr_percentile)
    rr_percentile = 50;
end

% Sample duration constants
sample_10ms  = round(fs * 0.01);
sample_70ms  = round(fs * 0.07);
sample_100ms = round(fs * 0.1);
sample_250ms = round(fs * 0.25);
sample_350ms = round(fs * 0.35);

N = length(ecg_rpeaks_index);
data_len = length(data_bias);
ecg_qrson_index  = ecg_qrson_index(:);
ecg_qrsoff_index = ecg_qrsoff_index(:);
ecg_rpeaks_index = ecg_rpeaks_index(:);

% === Prepare T-wave signal (same as fiducial_det_lsim lines 1146-1161) ===
ecg_Twave = data_bias;
new_baseline = nan(N - 2, 1);
for p = 2:N-1
    idx_qrs = ecg_qrson_index(p):ecg_qrsoff_index(p);
    if isempty(idx_qrs) || idx_qrs(1) < 1 || idx_qrs(end) > data_len
        continue;
    end
    ecg_Twave(idx_qrs) = linspace(ecg_Twave(idx_qrs(1)), ...
        ecg_Twave(idx_qrs(end)), length(idx_qrs));
    new_baseline(p - 1) = ecg_Twave(idx_qrs(end));
end
ecg_Twave = lp_filter_zero_phase(ecg_Twave, 10 / fs);
ecg_Twave = ecg_Twave - median(new_baseline, 'omitnan');

% === Process each cluster ===
% Pre-allocate so every cluster index is always valid downstream,
% even when a cluster has too few beats to estimate QT/RR.
for c = length(num_cls):-1:1
    qt_rr_prior(c) = empty_cluster_struct(num_cls(c));
end

for c = 1:length(num_cls)

    % Inner-beat indices for this cluster (1-based for beats 2..N-1)
    ind_c = index_clustering(index_clustering(:, 2) == num_cls(c), 1);
    p_all = ind_c + 1;  % corresponding ecg_rpeaks_index positions

    % Quality and avg RR for beats in this cluster
    quality_c = beat_quality_twave(p_all);
    avg_rr_c  = avg_intervals_ecg(p_all);

    % --- Select clean subset ---
    % Step 1: high quality beats
    q_mask = quality_c >= max(quality_thr, prctile(quality_c, 90));

    % Step 2: among high quality, prefer lower heart rate (higher avg_intervals_ecg)
    if sum(q_mask) > 3
        rr_thr = prctile(avg_rr_c(q_mask), rr_percentile);
        clean_mask = q_mask & (avg_rr_c >= rr_thr);
    else
        clean_mask = q_mask;
    end

    % Ensure minimum beats per cluster (all clusters must be represented)
    min_beats = max(3, round(0.1 * length(ind_c)));
    if sum(clean_mask) < min_beats
        clean_mask = q_mask;
    end
    if sum(clean_mask) < min_beats
        % Relax quality threshold
        clean_mask = quality_c >= max(0.5, quality_thr - 0.25);
    end
    if sum(clean_mask) < 3
        % Fallback: use all beats in cluster
        clean_mask = true(size(ind_c));
    end

    clean_p = p_all(clean_mask);
    n_clean = length(clean_p);

    % --- Extract T-wave segments for clean beats ---
    twave_beats   = nan(n_clean, 2 * sample_350ms);
    twave_lengths = zeros(n_clean, 1);
    qrs_dur       = zeros(n_clean, 1);

    for k = 1:n_clean
        p = clean_p(k);
        if p + 1 > N, continue; end

        qrs_dur(k) = ecg_qrsoff_index(p) - ecg_qrson_index(p);

        t_end = ecg_qrsoff_index(p) + min( ...
            (sample_350ms + sample_70ms) * max(1, avg_intervals_ecg(p) / fs), ...
            ecg_qrson_index(p + 1) - ecg_qrsoff_index(p) - sample_100ms);
        seg_idx = ecg_qrsoff_index(p) + 1 : t_end;

        if isempty(seg_idx) || seg_idx(end) > data_len || seg_idx(1) < 1
            continue;
        end
        if length(seg_idx) < 2 * sample_70ms
            seg_idx = ecg_qrsoff_index(p) + 1 : min(data_len, ecg_qrsoff_index(p) + sample_250ms);
        elseif length(seg_idx) > 2 * sample_350ms
            seg_idx = seg_idx(1 : 2 * sample_350ms);
        end
        if seg_idx(end) > data_len, continue; end

        seg = ecg_Twave(seg_idx);
        twave_beats(k, 1:length(seg)) = seg;
        twave_lengths(k) = length(seg);
    end

    valid_beat = twave_lengths > 0;
    n_valid = sum(valid_beat);

    if n_valid < 2
        qt_rr_prior(c) = empty_cluster_struct(num_cls(c));
        continue;
    end

    % --- Build average T-wave template ---
    avg_len = round(mean(twave_lengths(valid_beat)));
    avg_len = max(avg_len, 2 * sample_70ms);
    avg_tw  = mean(twave_beats(valid_beat, 1:avg_len), 1, 'omitmissing');
    avg_tw  = movmean(avg_tw, [sample_10ms, sample_10ms]);
    avg_tw  = fillmissing(avg_tw, "nearest");

    % --- Detect T-peak on template (same logic as fiducial_det_lsim) ---
    skip  = 3 * sample_10ms;
    tsrch = avg_tw(skip + 1 : end);

    [~, Pv_max] = islocalmax(tsrch, 'MaxNumExtrema', 2, ...
        'MinSeparation', sample_70ms, 'MinProminence', 0.05 * std(tsrch));
    [~, Pv_min] = islocalmin(tsrch, 'MaxNumExtrema', 2, ...
        'MinSeparation', sample_70ms, 'MinProminence', 0.05 * std(tsrch));

    [mx_v, mx_i] = max(Pv_max);
    [mn_v, mn_i] = max(Pv_min);

    if mx_v > 0 && mn_v > 0
        if mx_v >= mn_v
            t_peak_tmpl = mx_i + skip;
        else
            t_peak_tmpl = mn_i + skip;
        end
    elseif mx_v > 0
        t_peak_tmpl = mx_i + skip;
    elseif mn_v > 0
        t_peak_tmpl = mn_i + skip;
    else
        bl = mean([avg_tw(1:sample_10ms), avg_tw(max(1,end-sample_10ms+1):end)]);
        [~, t_peak_tmpl] = max(abs(avg_tw - bl));
    end

    % --- Estimate Toff on template ---
    toff_tmpl = detect_toff(avg_tw, t_peak_tmpl, ...
        sample_10ms, sample_70ms, win_sample_T);

    qrs_dur_med = round(median(qrs_dur(valid_beat), 'omitnan'));
    qt_template = qrs_dur_med + toff_tmpl;

    % --- Per-beat QT and RR estimation ---
    qt_vals = nan(n_clean, 1);
    rr_vals = nan(n_clean, 1);

    for k = 1:n_clean
        if ~valid_beat(k), continue; end
        p = clean_p(k);
        L = twave_lengths(k);

        bsig = twave_beats(k, 1:L);
        bsig = movmean(bsig, [sample_10ms, sample_10ms]);
        bsig = fillmissing(bsig, "nearest");

        % Refine T-peak near template estimate
        t_pk = min(t_peak_tmpl, L - sample_10ms);
        rgn  = max(1, t_pk - 3*sample_10ms) : min(L, t_pk + 3*sample_10ms);
        bl_b = mean([bsig(1:min(L,sample_10ms)), bsig(max(1,L-sample_10ms+1):L)]);
        [~, loc] = max(abs(bsig(rgn) - bl_b));
        t_pk = rgn(loc);

        % Detect Toff on individual beat
        toff_beat = detect_toff(bsig, t_pk, ...
            sample_10ms, sample_70ms, win_sample_T);

        qt_vals(k) = qrs_dur(k) + toff_beat;
        rr_vals(k) = rr_intervals_ecg(min(p, length(rr_intervals_ecg)));
    end

    % Reject outlier QT values (robust MAD-based)
    ok = ~isnan(qt_vals);
    if sum(ok) > 5
        qt_md  = median(qt_vals(ok));
        qt_mad = 1.4826 * median(abs(qt_vals(ok) - qt_md));
        qt_vals(abs(qt_vals - qt_md) > 3 * qt_mad) = nan;
    end

    final = ~isnan(qt_vals) & ~isnan(rr_vals);

    qt_rr_prior(c).QT             = 1000*qt_vals(final)/fs;
    qt_rr_prior(c).RR             = 1000*rr_vals(final)/fs;
    qt_rr_prior(c).median_QT      = median(1000*qt_vals(final)/fs, 'omitnan');
    qt_rr_prior(c).median_RR      = median(1000*rr_vals(final)/fs, 'omitnan');
    qt_rr_prior(c).qt_template    = 1000*qt_template/fs;
    qt_rr_prior(c).num_beats      = sum(final);
    qt_rr_prior(c).cluster_id     = num_cls(c);
    qt_rr_prior(c).clean_beat_idx = clean_p(final);
end

% Clusters with too few beats inherit prior values from the best cluster,
% so downstream code can index qt_rr_prior(c) for every cluster c.
num_beats_arr = [qt_rr_prior.num_beats];
[max_beats, best_c] = max(num_beats_arr);
if max_beats >= 2
    for c = 1:length(num_cls)
        if qt_rr_prior(c).num_beats < 2 || isnan(qt_rr_prior(c).median_QT)
            cid = qt_rr_prior(c).cluster_id;
            qt_rr_prior(c) = qt_rr_prior(best_c);
            qt_rr_prior(c).cluster_id = cid;
        end
    end
end

end


%% ========================================================================
function toff = detect_toff(sig, t_peak, sample_10ms, sample_70ms, win_T)
%DETECT_TOFF  Detect T-wave offset using derivative-threshold after T-peak
%   Uses the same feature types as fiducial_det_lsim (abs-diff + moving std)
%   with Gaussian weighting to suppress late activity (P-wave of next beat).

    L = length(sig);
    after_start = min(t_peak + win_T, L);

    if after_start >= L - 2 * sample_10ms
        toff = L;
        return;
    end

    after = sig(after_start:end);
    n = length(after);

    if n < 3 * sample_10ms
        toff = L;
        return;
    end

    % Gaussian weighting: descending from peak (same convention as main code)
    gw_full = gausswin(4 * max(sample_70ms, n) + 1)';
    mid = ceil(length(gw_full) / 2);
    gw = gw_full(mid : min(length(gw_full), mid + n - 1));
    if length(gw) < n
        gw = [gw, gw(end) * ones(1, n - length(gw))];
    end

    % Derivative feature (absolute difference)
    nd = zeros(1, n);
    valid_rng = win_T + 1 : n;
    if ~isempty(valid_rng)
        nd(valid_rng) = abs(after(valid_rng) - after(valid_rng - win_T));
    end

    % Moving std feature
    ls = movstd(after, [win_T, win_T]);

    % Apply Gaussian weighting
    nd = nd .* gw;
    ls = ls .* gw;

    % Combined feature
    feat = (nd + ls) / 2;

    % Find where feature drops below threshold
    search_from = min(n, 2 * sample_10ms + 1);
    peak_f = max(feat(1 : min(n, 5 * sample_10ms)));

    if peak_f < 1e-10
        toff = after_start + round(n / 3);
        toff = min(toff, L);
        return;
    end

    thr = 0.20 * peak_f;
    idx = find(feat(search_from:end) < thr, 1, 'first');

    if isempty(idx)
        [~, idx] = min(feat(search_from:end));
    end
    idx = idx + search_from - 1;

    toff = after_start + idx - 1;
    toff = min(toff, L);
end


%% ========================================================================
function s = empty_cluster_struct(cid)
%EMPTY_CLUSTER_STRUCT  Default struct when cluster has too few valid beats
    s.QT             = [];
    s.RR             = [];
    s.median_QT      = nan;
    s.median_RR      = nan;
    s.qt_template    = nan;
    s.num_beats      = 0;
    s.cluster_id     = cid;
    s.clean_beat_idx = [];
end
