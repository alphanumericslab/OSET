% function test_ecg_atrial_ventricular_bss.m
clear; close all; clc;

datafilepath = '../../../DataFiles/physionet.org/files/ptbdb/1.0.0/'; fs = 1000.0;
% datafilepath = '../../../DataFiles/physionet.org/files/qtdb/1.0.0/'; fs = 250.0;

filelist = dir(fullfile([datafilepath, '**/*.mat']));  % get list of all mat files
fl = 0.1;
fh = 80.0;
f_hr = 1.4; % default initial heart rate value in Hz
top_source_num = 4;
num_rounds = 2; % default number of R-peak detection rounds; see peak_det_local_search() help for details
baseline_fraction = 0.1; % percentage length of the beat (from the start) considered as the baseline segment (roughly)
avg_beat_right_expansion_frac = 1.2;
% refch = 1;
% data_augmented = [];
for k = 54 : 54%length(filelist)
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val;
    headerfilename = [filelist(k).folder '/' filelist(k).name];
    header_ID = fopen([headerfilename(1:end-3) 'hea'], 'r');
    line_ = split(fgetl(header_ID), ' ');
    num_ch = str2double(line_{2});
    ch_names = cell(1, num_ch);
    for ch = 1 : num_ch
        line_ = split(fgetl(header_ID), ' ');
        ch_names{ch} = line_{end};
    end

    data = data - lp_filter_zero_phase(data, fl/fs);
    data = lp_filter_zero_phase(data, fh/fs);

    % bl_params.wlen1 = round(fs*0.7);
    % bl_params.wlen2 = round(fs*0.65);
    % bl_params.wlen3 = round(fs*0.75);
    % bl_params.wlen4 = round(fs*0.7);
    % baseline = baseline_filter(data, 'MDMDMDMN', bl_params);
    % data = data - baseline;

    data = data(:, round(5*fs) : round(20*fs));

    P_wave_indexes = [];
    T_wave_indexes = [];
    position = cell(size(data, 1),1);
    peak_indexes = cell(size(data, 1),1);
    for ch = 1 : size(data, 1)
        %{
        peak_detector_params.saturate = 1;
        peak_detector_params.k_sigma = 4;
        peak_detector_params.hist_search_th = 0.9;
        peak_detector_params.rpeak_search_wlen = 0.4; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
        peak_detector_params.filter_type = 'MDMN';%'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER', 'MDMN', 'WAVELET'
        [~, peak_indexes, ~] = peak_det_probabilistic(data(ch, :), fs, peak_detector_params);
        %}

        [~, peak_indexes{ch}] = peak_det_local_search(data(ch, :), f_hr/fs, [], num_rounds);

        data_ch = data(ch, :);
        heasig.nsig = 1; % N is number of channels
        heasig.freq = fs;
        heasig.nsamp = length(data_ch);
        [position{ch}, ~, ~] = wavedet_3D(data_ch(:), peak_indexes{ch}, heasig);

        for p = 1 : length(position{ch}.Pon)
            if ~isnan(position{ch}.Pon(p)) && ~isnan(position{ch}.Poff(p))
                P_wave_indexes = cat(2, P_wave_indexes, position{ch}.Pon(p) : position{ch}.Poff(p));
                % P_wave_indexes = cat(2, P_wave_indexes, position{ch}.Pon(p) : position{ch}.QRSon(p)-1);
                % P_wave_indexes = cat(2, P_wave_indexes, position{ch}.Pon(p) : position{ch}.qrs(p)-1);
            end
        end

        for t = 1 : length(position{ch}.Ton)
            if ~isnan(position{ch}.Ton(t)) && ~isnan(position{ch}.Toff(t))
                T_wave_indexes = cat(2, T_wave_indexes, position{ch}.Ton(t) : position{ch}.Toff(t));
                % T_wave_indexes = cat(2, T_wave_indexes, position{ch}.qrs(t) : position{ch}.Toff(t));
            end
        end
    end

    P_wave_indexes = unique(P_wave_indexes);
    T_wave_indexes = unique(T_wave_indexes);

    % [s, W, A, B, C] = nonstationary_component_analysis(data, P_wave_indexes, 1 : length(data_ch));

    [s_nsca, W, A, B, C, lambda] = nonstationary_component_analysis(data, P_wave_indexes, T_wave_indexes);
    data_atrial = A(:, 1 : top_source_num)*s_nsca(1 : top_source_num, :);
    data_ventricular = A(:, end - top_source_num + 1: end)*s_nsca(end - top_source_num + 1: end, :);
    %{
    % Equivalent as above using MICA
    [Proj_tilde, ~] = mica_projectors(A);
    Map_atrial = zeros(size(data, 1));
    Map_ventricular = zeros(size(data, 1));
    for mm = 1 : top_source_num
        Map_atrial = Map_atrial + Proj_tilde{mm};
        Map_ventricular = Map_ventricular + Proj_tilde{end-mm+1};
    end
    data_atrial = Map_atrial * data;
    data_ventricular = Map_ventricular * data;
    %}

    % % [s_atrial, W_atrial, A_atrial, B_atrial, C_atrial] = nonstationary_component_analysis(data, P_wave_indexes, 1 : length(data_ch));
    % [s_atrial, W_atrial, A_atrial, B_atrial, C_atrial] = nonstationary_component_analysis(data, P_wave_indexes, setdiff(1 : length(data_ch), P_wave_indexes));
    % data_atrial = A_atrial(:, 1 : top_source_num)*s_atrial(1 : top_source_num, :);
    %
    % % [s_ventricular, W_ventricular, A_ventricular, B_ventricular, C_ventricular] = nonstationary_component_analysis(data, T_wave_indexes, 1 : length(data_ch));
    % [s_ventricular, W_ventricular, A_ventricular, B_ventricular, C_ventricular] = nonstationary_component_analysis(data, T_wave_indexes, setdiff(1 : length(data_ch), T_wave_indexes));
    % data_ventricular = A_ventricular(:, 1 : top_source_num)*s_ventricular(1 : top_source_num, :);

    % W = jadeR(data - mean(data, 2));
    % data_jade = W * data;

    % data_jade = data_jade - mean(data_jade, 2);
    % data_jade = normr(data_jade);

    plot_multichannel_data(data, 5, 'k', fs)
    plot_multichannel_data(s_nsca, 5, 'b', fs)
    
    % plot_multichannel_data(s_atrial, 5, 'b', fs)
    % plot_multichannel_data(s_ventricular, 5, 'r', fs)

    % plot sources in the source space
    refch = 1;
    beat_width = round(1.0 * median(diff([1, peak_indexes{refch}, size(data, 2)]))); % Add the first and last indexes
    if mod(beat_width, 2) == 0
        beat_width = beat_width + 1;
    end
    beat_left_wing = round(1.0 * (beat_width - 1) / 2);
    beat_right_wing = round(avg_beat_right_expansion_frac*(beat_width - 1) / 2);
    beat_width_expanded = beat_left_wing + beat_right_wing + 1;
    ECG_robust_median_source_space = zeros(size(s_nsca, 1), beat_width_expanded);
    for source = 1 : size(s_nsca, 1)
        stacked_beats = event_stacker(s_nsca(source, :), peak_indexes{refch}, [beat_left_wing, beat_right_wing]); % stack the beats
        dc_level = mean(mean(stacked_beats(:, 1:round(baseline_fraction * beat_width_expanded)))); % estimate the DC level from the baseline segment
        stacked_beats = stacked_beats - dc_level; % remove the baseline level offset
        [~, ~, ECG_robust_median_source_space(source, :)] = robust_weighted_average(stacked_beats); % find the robust weighted median beat
        % tm = (0 : size(stacked_beats, 2) - 1) / fs;
        % figure
        % plot(tm, stacked_beats', 'color', 0.7*ones(1, 3))
        % hold on
        % plot(tm, ECG_robust_median_source_space, 'b', 'linewidth', 2)
        % title(['Source #' num2str(source)]);
        % xlabel('time(s)');
        % ylabel('Amplitude(normalized)');
        % grid
        % set(gca, 'fontsize', 18)
    end
    tm = (0 : beat_width_expanded - 1) / fs;
    figure
    plot(tm, ECG_robust_median_source_space')
    title('Source space average beats');
    legend();
    xlabel('time(s)');
    ylabel('Amplitude(normalized)');
    grid
    set(gca, 'fontsize', 18)
    legend(string(1:size(s_nsca, 1)))


    % [~, peak_indexes] = peak_det_local_search(data(refch, :), f_hr/fs, [], num_rounds);

    P_onset = zeros(size(data, 1), 1);
    P_offset = zeros(size(data, 1), 1);
    QRS_onset = zeros(size(data, 1), 1);
    QRS_offset = zeros(size(data, 1), 1);
    T_onset = zeros(size(data, 1), 1);
    T_offset = zeros(size(data, 1), 1);
    for ch = 1 : size(data, 1)
        beat_width = round(1.0 * median(diff([1, peak_indexes{ch}, size(data, 2)]))); % Add the first and last indexes
        if mod(beat_width, 2) == 0
            beat_width = beat_width + 1;
        end
        beat_left_wing = round(1.0 * (beat_width - 1) / 2);
        beat_right_wing = round(avg_beat_right_expansion_frac*(beat_width - 1) / 2);
        beat_width_expanded = beat_left_wing + beat_right_wing + 1;

        P_onset(ch) = round(beat_left_wing + 1 - median(position{ch}.qrs - position{ch}.Pon, "omitmissing"));
        P_offset(ch) = round(beat_left_wing + 1 - median(position{ch}.qrs - position{ch}.Poff, "omitmissing"));
        QRS_onset(ch) = round(beat_left_wing + 1 - median(position{ch}.qrs - position{ch}.QRSon, "omitmissing"));
        QRS_offset(ch) = round(beat_left_wing + 1 - median(position{ch}.qrs - position{ch}.QRSoff, "omitmissing"));
        T_onset(ch) = round(beat_left_wing + 1 - median(position{ch}.qrs - position{ch}.Ton, "omitmissing"));
        T_offset(ch) = round(beat_left_wing + 1 - median(position{ch}.qrs - position{ch}.Toff, "omitmissing"));

        [~, peak_indexes_atrial] = peak_det_local_search(data_atrial(ch, :), f_hr/fs, [], num_rounds);
        [~, peak_indexes_ventricular] = peak_det_local_search(data_ventricular(ch, :), f_hr/fs, [], num_rounds);

        hr = 60 * fs ./ diff(peak_indexes{ch});
        hr_atrial = 60 * fs ./ diff(peak_indexes_atrial);
        hr_ventricular = 60 * fs ./ diff(peak_indexes_ventricular);

        % calculate the time domain average beat
        stacked_beats = event_stacker(data(ch, :), peak_indexes{ch}, [beat_left_wing, beat_right_wing]); % stack the beats
        dc_level = mean(mean(stacked_beats(:, 1:round(baseline_fraction * beat_width_expanded)))); % estimate the DC level from the baseline segment
        stacked_beats = stacked_beats - dc_level; % remove the baseline level offset
        [~, ~, ECG_robust_median] = robust_weighted_average(stacked_beats); % find the robust weighted median beat

        % % calculate the time domain average beat
        % stacked_beats_jade = event_stacker(data_jade(ch, :), peak_indexes, [beat_left_wing, beat_right_wing]); % stack the beats
        % dc_level_jade = 0;   %     mean(mean(stacked_beats_jade(:, 1:round(baseline_fraction * beat_width_expanded)))); % estimate the DC level from the baseline segment
        % stacked_beats_jade = stacked_beats_jade - dc_level_jade; % remove the baseline level offset
        % [~, ~, ECG_robust_median_jade] = robust_weighted_average(stacked_beats_jade); % find the robust weighted median beat

        % calculate the time domain average beat
        stacked_beats_atrial = event_stacker(data_atrial(ch, :), peak_indexes{ch}, [beat_left_wing, beat_right_wing]); % stack the beats
        dc_level_atrial =  mean(mean(stacked_beats_atrial(:, 1:round(baseline_fraction * beat_width_expanded)))); % estimate the DC level from the baseline segment
        stacked_beats_atrial = stacked_beats_atrial - dc_level_atrial; % remove the baseline level offset
        [~, ~, ECG_robust_median_atrial] = robust_weighted_average(stacked_beats_atrial); % find the robust weighted median beat

        % calculate the time domain average beat
        stacked_beats_ventricular = event_stacker(data_ventricular(ch, :), peak_indexes{ch}, [beat_left_wing, beat_right_wing]); % stack the beats
        dc_level_ventricular =  mean(mean(stacked_beats_ventricular(:, 1:round(baseline_fraction * beat_width_expanded)))); % estimate the DC level from the baseline segment
        stacked_beats_ventricular = stacked_beats_ventricular - dc_level_ventricular; % remove the baseline level offset
        [~, ~, ECG_robust_median_ventricular] = robust_weighted_average(stacked_beats_ventricular); % find the robust weighted median beat

        time = (0 : size(data, 2) - 1) / fs;
        lgnd = {};
        figure
        subplot(311)
        plot(time(peak_indexes{ch}(2:end)), hr); lgnd = cat(1, lgnd, 'HR');
        hold on
        plot(time(peak_indexes_atrial(2:end)), hr_atrial); lgnd = cat(1, lgnd, 'atrial-basd HR');
        plot(time(peak_indexes_ventricular(2:end)), hr_ventricular); lgnd = cat(1, lgnd, 'ventricular-based HR');
        legend(lgnd);
        xlabel('time(s)');
        ylabel('HR (BPM)');
        grid
        set(gca, 'fontsize', 18)

        lgnd = {};
        subplot(312)
        plot(time, data(ch, :)); lgnd = cat(1, lgnd, 'ECG');
        hold on
        plot(time(peak_indexes{ch}), data(ch, peak_indexes{ch}), 'ro', 'markersize', 16); lgnd = cat(1, lgnd, 'R-peaks');
        try plot(time(position{ch}.Pon), data(ch, position{ch}.Pon), 'k>', 'markersize', 16); lgnd = cat(1, lgnd, 'P-onset');
        catch
        end
        try plot(time(position{ch}.Poff), data(ch, position{ch}.Poff), 'k<', 'markersize', 16); lgnd = cat(1, lgnd, 'P-offset');
        catch
        end
        try plot(time(position{ch}.QRSon), data(ch, position{ch}.QRSon), 'k>', 'markersize', 16); lgnd = cat(1, lgnd, 'QRS-onset');
        catch
        end
        try plot(time(position{ch}.QRSoff), data(ch, position{ch}.QRSoff), 'k<', 'markersize', 16); lgnd = cat(1, lgnd, 'QRS-offset');
        catch
        end
        try plot(time(position{ch}.Ton), data(ch, position{ch}.Ton), 'k>', 'markersize', 16); lgnd = cat(1, lgnd, 'T-onset');
        catch
        end
        try plot(time(position{ch}.Toff), data(ch, position{ch}.Toff), 'k<', 'markersize', 16); lgnd = cat(1, lgnd, 'T-offset');
        catch
        end
        plot(time, data_atrial(ch, :));  lgnd = cat(1, lgnd, 'atrial');
        plot(time, data_ventricular(ch, :));  lgnd = cat(1, lgnd, 'ventricular');
        legend(lgnd);
        xlabel('time(s)');
        ylabel('Amplitude(uV)');
        grid
        set(gca, 'fontsize', 18)

        tm = (0 : size(stacked_beats, 2) - 1) / fs;
        subplot(325)
        plot(tm, stacked_beats', 'color', 0.7*ones(1, 3))
        hold on
        plot(tm, ECG_robust_median, 'k', 'linewidth', 2)
        try plot(tm(P_onset(ch)), ECG_robust_median(P_onset(ch)), 'r>', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(P_offset(ch)), ECG_robust_median(P_offset(ch)), 'r<', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(beat_left_wing + 1), ECG_robust_median(beat_left_wing + 1), 'bo', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(QRS_onset(ch)), ECG_robust_median(QRS_onset(ch)), 'g>', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(QRS_offset(ch)), ECG_robust_median(QRS_offset(ch)), 'g<', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(T_onset(ch)), ECG_robust_median(T_onset(ch)), 'b>', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(T_offset(ch)), ECG_robust_median(T_offset(ch)), 'b<', 'linewidth', 1, 'markersize', 12);
        catch
        end
        xlabel('time(s)');
        ylabel('Amplitude(uV)');
        title('overall median beat')
        grid
        set(gca, 'fontsize', 18)

        tm = (0 : size(stacked_beats_atrial, 2) - 1) / fs;
        subplot(326)
        plot(tm, stacked_beats_atrial', 'color', [0.6, 0.6, 0.9])
        hold on
        % xlabel('time(s)');
        % ylabel('Amplitude(uV)');
        % title('atrial median beat')
        % grid
        % set(gca, 'fontsize', 18)

        %tm = (0 : size(stacked_beats_ventricular, 2) - 1) / fs;
        % subplot(236)
        plot(tm, stacked_beats_ventricular', 'color', [0.9, 0.6, 0.6])
        % plot(tm, norm(stacked_beats_ventricular(:)) * stacked_beats_jade' / norm(stacked_beats_jade(:)), 'color', [0.6, 0.9, 0.6])
        % hold on
        h_av = plot(tm, ECG_robust_median_atrial + ECG_robust_median_ventricular, 'k', 'linewidth', 3);
        h_a = plot(tm, ECG_robust_median_atrial, 'b', 'linewidth', 2);
        h_v = plot(tm, ECG_robust_median_ventricular, 'r', 'linewidth', 2);
        try plot(tm(P_onset(ch)), 0, 'k>', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(P_offset(ch)), 0, 'k<', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(beat_left_wing + 1), ECG_robust_median_atrial(beat_left_wing + 1), 'ko', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(beat_left_wing + 1), ECG_robust_median_ventricular(beat_left_wing + 1), 'ko', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(QRS_onset(ch)), 0, 'k>', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(QRS_offset(ch)), 0, 'k<', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(T_onset(ch)), 0, 'k>', 'linewidth', 1, 'markersize', 12);
        catch
        end
        try plot(tm(T_offset(ch)), 0, 'k<', 'linewidth', 1, 'markersize', 12);
        catch
        end
        % h_jade = plot(tm, norm(ECG_robust_median_ventricular) * ECG_robust_median_jade / norm(ECG_robust_median_jade), 'g', 'linewidth', 2);
        % legend([h_a, h_v, h_jade], {'atrial', 'ventricular', 'jade'});
        legend([h_a, h_v, h_av], {'atrial', 'ventricular', 'atrial+ventricular'});
        xlabel('time(s)');
        ylabel('Amplitude(uV)');
        title('atrial and ventricular median beat')
        grid
        set(gca, 'fontsize', 18)

        sgtitle(['Record ', filelist(k).name ', Channel ', ch_names{ch}], 'interpreter', 'none', 'fontsize', 24)
    end
end
