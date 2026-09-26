# LSIM-Deli: ECG delineation with Latent Structure Influence Models

LSIM-Deli finds the fiducial points of every beat in a single-lead ECG: the onset, peak, and offset of the P-wave, QRS complex, and T-wave. There are two delineators, with the same syntax and outputs:

- [`ecg_delineate_lsim.m`](ecg_delineate_lsim.m) — the main delineator. Use it for routine, short-to-moderate recordings (resting ECG, a few minutes) where the heart rate stays fairly stable.
- [`ecg_delineate_lsim_long.m`](ecg_delineate_lsim_long.m) — same algorithm, plus a per-cluster QT/RR prior that rescales the T-wave search window to each beat's own RR interval. Use it for **Holter recordings and exercise/stress-test ECG**, where the heart rate drifts substantially over the recording and a fixed search window makes Toff drift with it.

Both share every helper function except the T-wave search window, so switching between them only means renaming the call.

| Wave | Fiducial points |
|------|-----------------|
| P    | `Pon`, `P`, `Poff` |
| QRS  | `QRSon`, `R`, `QRSoff` |
| T    | `Ton`, `T`, `Toff` |

## How it works

1. **Beat preparation** (`build_beat_context`). The signal is high-passed at 0.5 Hz and the R-peaks are cleaned: a beat is dropped if its correlation-based quality score is low or its RR interval is not physiological. The remaining beats are grouped into morphology clusters with fuzzy c-means. The number of clusters is chosen by the Xie–Beni index.
2. **QRS** (`detect_qrs`). Each beat has one segment before the R-peak and one after it. From each segment, two features are taken: the absolute local difference and the moving standard deviation. An LSIM model decodes where the signal switches between active and isoelectric states, which gives `QRSon` and `QRSoff`.
3. **T-wave** (`detect_twave`). The QRS is cut out, the signal is low-passed at 10 Hz, and an average T-wave is built for each cluster. That template sets the T-wave polarity (positive, negative, or biphasic) and a starting estimate of the peak. `Ton` and `Toff` then come from LSIM on/off decoding of each beat. In `ecg_delineate_lsim_long`, the segment searched for Toff is first rescaled per beat from a QT/RR prior (`compute_qt_rr_prior`, see below) instead of a fixed, heart-rate-independent window.
4. **P-wave** (`detect_pwave`). The P-wave is found the same way in the PR segment. Each beat also gets a `P_score` between 0 and 1, which can be used to drop unreliable P-waves.
5. **Post-processing** (`regularize_interval`, optional). Priors built from the RR intervals regularise the intervals within each cluster and correct outliers.

All three waves use the same LSIM transition decoder, `lsim_onoff_detection`, and the same two features.

`compute_qt_rr_prior` (used only by `ecg_delineate_lsim_long`) picks, per morphology cluster, the cleanest beats with the lowest heart rate, builds a T-wave template from them, and fits a QT/RR curve. Each beat's T-wave search window is then this curve evaluated at that beat's own smoothed RR interval, instead of a window sized from a fixed multiple of `sample_350ms`. This keeps the search window realistic across a wide range of heart rates within the same recording — the situation a short resting-ECG recording doesn't have, but a multi-hour Holter or an exercise/stress test does.

## Requirements

- MATLAB R2021a or newer. It needs the Signal Processing, Statistics and Machine Learning, and Fuzzy Logic toolboxes (FCM clustering). The Parallel Computing Toolbox is optional; `parfor` is used when it is available.
- OSET functions on the path: [`matlab/`](../../matlab), which includes `peak_det_likelihood_long_recs` and `lp_filter_zero_phase`.
- The LSIM toolbox (`em_lsim`), which ships with OSET in [`external/chmm-lsim-matlab-toolbox`](../../external/chmm-lsim-matlab-toolbox) and is also at <https://github.com/sajjadkarimi91/chmm-lsim-karimi-toolbox>.
- This folder, which contains helpers such as `em_interval_calc`, `lsim_fit_detect`, and (for `ecg_delineate_lsim_long`) `compute_qt_rr_prior`.

```matlab
oset_path = 'path/to/OSET';
addpath(genpath(fullfile(oset_path, 'matlab')))
addpath(genpath(fullfile(oset_path, 'external', 'chmm-lsim-matlab-toolbox')))
addpath(fullfile(oset_path, 'under-development', 'lsim-deli'))
```

> Avoid `addpath(genpath(oset_path))` on the whole toolbox. It also adds the `depricated/` folders, whose older copies can shadow the current functions.

## Quick start

```matlab
% ecg: 1 x T single-lead ECG, fs: sampling rate (Hz)

% Simplest call: R-peaks are detected internally
positions = ecg_delineate_lsim(ecg, fs);

% With your own R-peaks and post-processing on
[positions, EXITFLAG] = ecg_delineate_lsim(ecg, fs, rpeaks, 1);

% Internal R-peaks plus optional arguments
[positions, EXITFLAG] = ecg_delineate_lsim(ecg, fs, [], 1);

qt_ms = 1000 * (positions.Toff - positions.QRSon) / fs;   % per-beat QT interval
```

If you leave out the R-peaks, or pass `[]`, the function detects them itself with `peak_det_likelihood_long_recs`. It tries three detection thresholds and keeps the one with the lowest short-term HRV. The high-pass cutoff it used is returned in `positions.rpeak_bp_lower_cutoff_hz`.

[`sample_deli.m`](sample_deli.m) is a complete demo. It runs on every `.mat` record in this folder (for example `sel32.mat` from the QT database), plots the detected points over the expert annotations, and prints the mean absolute error for each fiducial point. It calls `ecg_delineate_lsim`; swap in `ecg_delineate_lsim_long` for Holter or exercise/stress-test data.

## Syntax

`ecg_delineate_lsim_long` has the identical syntax; replace the function name below.

```matlab
[positions, EXITFLAG] = ecg_delineate_lsim(data, fs)
[positions, EXITFLAG] = ecg_delineate_lsim(data, fs, ecg_rpeaks_index)
[positions, EXITFLAG] = ecg_delineate_lsim(data, fs, ecg_rpeaks_index, ...
    flag_post_processing, flag_prune_P, win_qrs, win_T, win_P, ...
    max_clusters, twave_shape, time_prior_mode)
```

Only `data` and `fs` are required. Pass `[]` for any optional input to keep its default. To use optional inputs 4–11 without your own R-peaks, pass `[]` as `ecg_rpeaks_index`.

| # | Input | Default | Description |
|---|-------|---------|-------------|
| 1 | `data` | required | Single-lead ECG vector. Remove powerline noise and baseline wander before calling. |
| 2 | `fs` | required | Sampling rate in Hz. |
| 3 | `ecg_rpeaks_index` | `[]` | R-peak positions in samples. If empty or left out, R-peaks are detected internally. |
| 4 | `flag_post_processing` | `0` | `1` turns on the RR-prior refinement. It also removes low-quality beats more strictly (quality below 0.5 instead of 0.4) and scales the feature windows to the heart rate. Recommended. |
| 5 | `flag_prune_P` | `0` | `1` removes P-waves with a low `P_score`. |
| 6 | `win_qrs` | `0.01` s | Feature window for QRS onset and offset. |
| 7 | `win_T` | `0.02` s | Feature window for T-wave onset and offset. |
| 8 | `win_P` | `0.01` s | Feature window for P-wave onset and offset. |
| 9 | `max_clusters` | `max(3, min(6, ceil(T/(10*fs))))` | Maximum number of FCM morphology clusters. |
| 10 | `twave_shape` | `'none'` | T-wave shape to enforce: `'none'` detects it automatically; the others are `'bi-phasic'`, `'min'`, and `'max'`. |
| 11 | `time_prior_mode` | `'disabled'` | `'disabled'` uses a scalar prior of 1. `'normalized'` weights the difference feature by the local moving standard deviation divided by its 90th percentile. |

## Outputs

`positions` has one entry per input R-peak, in samples. Beats removed during preprocessing are `NaN`.

| Field | Description |
|-------|-------------|
| `Pon`, `P`, `Poff` | P-wave onset, peak, and offset |
| `QRSon`, `R`, `QRSoff` | QRS onset, R-peak (the input R-peaks), and QRS offset |
| `Ton`, `T`, `Toff` | T-wave onset, first T-wave peak, and offset |
| `P_score` | P-wave quality in [0, 1] |
| `beat_quality_score` | Correlation-based beat quality `q`; `0` means the beat was removed |
| `beat_snr` | Beat SNR in dB, computed as `10*log10(q^2/(1-q^2))`. `\|q\|` is clamped to [1e-4, 0.9999], so a removed beat gives -100 dB. |
| `index_clustering` | FCM morphology cluster of each beat. It is `NaN` for the first and last beats and for removed beats. Label numbers can change between runs; the fiducial points do not. |
| `rpeak_bp_lower_cutoff_hz` | High-pass cutoff used for internal R-peak detection. It is `NaN` when you supplied the R-peaks. |

`EXITFLAG.status` is `'succeeded'` or `'failed'`. If a stage throws an error, the function still returns every point it found up to that stage and fills the rest with `NaN`. The `MException` is stored in `EXITFLAG.message`.

`ecg_delineate_lsim_long` returns the same fields. Compared to `ecg_delineate_lsim`, only `Toff` (and, through its search window sitting after the previous beat's Toff, occasionally `Pon` by a sample or two) can differ; `Pon`/`P`/`Poff`/`QRSon`/`R`/`QRSoff`/`Ton`/`T` are unaffected.

## Files

| File | Purpose |
|------|---------|
| [`ecg_delineate_lsim.m`](ecg_delineate_lsim.m) | **Main ECG delineator** |
| [`ecg_delineate_lsim_long.m`](ecg_delineate_lsim_long.m) | **Delineator for Holter / exercise / stress-test ECG**, with a QT/RR-prior-adjusted Toff search window |
| [`compute_qt_rr_prior.m`](compute_qt_rr_prior.m) | Per-cluster QT/RR prior, used only by `ecg_delineate_lsim_long` |
| [`sample_deli.m`](sample_deli.m) | Demo and evaluation script |
| [`lsim_fit_detect.m`](lsim_fit_detect.m) | Fits the LSIM model and decodes transitions |
| [`em_interval_calc.m`](em_interval_calc.m) | RR-interval estimate |
| [`fiducial_det_ppg.m`](fiducial_det_ppg.m) | Equivalent LSIM delineator for PPG |
| [`fiducial_det_lsim.m`](fiducial_det_lsim.m) | **Deprecated.** Forwards to `ecg_delineate_lsim` |
| [`depricated/`](depricated) | Archived earlier implementations (`fiducial_det_lsim_v1.m`, `fiducial_det_lsim_v2.m`) |

## Migrating from `fiducial_det_lsim`

`fiducial_det_lsim` now forwards to `ecg_delineate_lsim`. The first call in a session shows the warning `OSET:fiducial_det_lsim:deprecated`. Existing code keeps working without changes, and the output struct now has three extra fields: `beat_snr`, `index_clustering`, and `rpeak_bp_lower_cutoff_hz`.

To migrate, rename the call and **swap the second and third inputs**. `ecg_delineate_lsim` takes `fs` second, so that every input after `fs` is optional. The optional inputs keep their positions (4–11).

```matlab
positions = fiducial_det_lsim(ecg, rpeaks, fs, 1);    % old
positions = ecg_delineate_lsim(ecg, fs, rpeaks, 1);   % new
```

If you pass the old order to `ecg_delineate_lsim` by mistake, it stops with the error `ecg_delineate_lsim:InputOrder` because `fs` is not a scalar.

To hide the warning: `warning('off', 'OSET:fiducial_det_lsim:deprecated')`.

## Authors

Sajjad Karimi and Reza Sameni, 2024.
Part of [OSET, the Open-Source Electrophysiological Toolbox](https://github.com/alphanumericslab/OSET).
