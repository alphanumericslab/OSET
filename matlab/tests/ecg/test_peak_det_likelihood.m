% Test script for R-peak detection with peak_det_likelihood

close all
clear
clc

% Change this path to where you have the .mat data files
% Note: if needed, convert WFDB dat files in bash using this command:
%   find ./*/ -type f -execdir wfdb2mat -r {} \;

datafilepath = '../../../../../DataFiles/Physionet.org/files/qtdb/1.0.0/';fs = 250.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
%%% datafilepath = '../../../../../DataFiles/physionet.org/files/ptb-xl/1.0.3/records500/'; fs = 500.0;
% datafilepath = '../../../../../DataFiles/physionet.org/files/ptb-xl/1.0.3/records100/'; fs = 100.0;
% datafilepath = '../../../../../DataFiles/physionet.org/files/ptbdb/1.0.0/'; fs = 1000.0;
% datafilepath = '../../../../../DataFiles/physionet.org/files/mitdb/1.0.0/'; fs = 360.0;

PLOT_FIGURES = false; % plot figures or not
SAVE_FIGURES = false; % save r-peak figures or not (figure is not displayed if this flag is active)
SAVE_RESULTS = true; % save r-peak results or not


output_results_folder = './results/';
if ~exist(output_results_folder, 'dir') && (SAVE_FIGURES || SAVE_RESULTS)
    % Folder does not exist so create it
    mkdir(output_results_folder);
end

% peak_det_likelihood has many internal parameters and flags, which all
% have default values. the following lines are only examples (they can be
% commented out, since they match the internal parameter defaults)
peak_detector_params.RETURN_SIGNAL_PEAKS = true; % return signal peaks or energy envelope peaks
peak_detector_params.PLOT_RESULTS = false; % plot the results using the internal plot function of peak_det_likelihood or not
peak_detector_params.PLOT_DIAGNOSTIC = false; % diagnostic mode (do not activate unless diving deep into the code! run only on short segments, since many figures are created)
peak_detector_params.verbose = false; % reports all the default values for the internal parameters of peak_det_likelihood, which can be modified through this data structure if needed.
%%%% peak_detector_params.REFINE_PEAKS = false;

overlap_time = 1.0; % overlap between segments for continuity (1.0-2.0 seconds is enough)
seg_len_time = 10.0; % segment length in seconds

filelist = dir(fullfile([datafilepath, '**/*.mat']));  % get list of all mat files of interest
for k = 1 : length(filelist) % Sweep over all or specific records
    % TSTART = tic;
    record_name = filelist(k).name;
    % if isequal(record_name, "sel14046m.mat")
    datafilename = fullfile(filelist(k).folder, record_name);
    full_fname = string(record_name(1:end-4));

    data = load(datafilename);
    data = data.val;

    sig_len = size(data, 2);

    header_ID = fopen([datafilename(1:end-3) 'hea'], 'r');
    line_ = split(fgetl(header_ID), ' ');
    num_ch = str2double(line_{2});
    ch_names = cell(1, num_ch);

    for ch = 1 : num_ch  % sweep over all or a single desired channel
        line_ = split(fgetl(header_ID), ' ');
        ch_names{ch} = line_{end}; %'visible','off'

        % 90.0s is no magic number below; for long signals it's better to use
        % peak_det_likelihood_long_recs, which internally breaks the
        % input, runs peak_det_likelihood and stiches the results
        % together
        if (sig_len/fs) < 90.0
            [peaks, peaks_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood(data(ch, :), fs, peak_detector_params);
        else
            [peaks, peaks_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(data(ch, :), fs, seg_len_time, overlap_time, peak_detector_params);
        end

        hr = fs * 60.0 ./ diff(peaks_indexes); % HR in bpm
        hr_consensus = fs * 60.0 ./ diff(peak_indexes_consensus); % HR in bpm

        if PLOT_FIGURES
            lgnds = {};
            tm = (0 : sig_len - 1) / fs;
            if SAVE_FIGURES
                fig = figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8], 'visible','off');
            else
                fig = figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);

            end
            subplot(3, 1, [1, 2])
            scatter(tm, data(ch, :), 75, qrs_likelihood'*[0.5, 0, 0] + 0.5, 'filled'); lgnds = cat(2, lgnds, 'QRS likelihood (color-coded from gray to red)');
            hold on
            plot(tm , data(ch, :), 'b', 'linewidth', 2); lgnds = cat(2, lgnds, 'signal');
            plot(tm(peaks_indexes), data(ch, peaks_indexes), 'ro', 'MarkerFaceColor','r', 'markersize', 12); lgnds = cat(2, lgnds, 'R-peaks detected');
            plot(tm(peak_indexes_consensus), data(ch, peak_indexes_consensus), 'go', 'MarkerFaceColor','g', 'markersize', 11);  lgnds = cat(2, lgnds, 'R-peaks corrected');
            grid
            set(gca, 'fontsize', 18, 'box', 'on')
            xlim([tm(1), tm(end)])
            xlabel('time[s]');
            ylabel('Amplitude');
            legend(lgnds, 'orientation', 'horizontal');
            title(['ECG and R-peaks of ', record_name, ', lead ', ch_names{ch}], 'interpreter', 'none');

            subplot(3, 1, 3)
            hold on
            plot(tm(peaks_indexes(2:end)), hr)
            plot(tm(peak_indexes_consensus(2:end)), hr_consensus)
            grid
            xlabel('time[s]');
            ylabel('beats per minute (bpm)');
            title('Heart rate');
            legend({'heart rate', 'heart rate corrected'}, 'orientation', 'horizontal')
            set(gca, 'fontsize', 18, 'box', 'on')
            % xlim([tm(peaks_indexes(2)), tm(peaks_indexes(end))])
            xlim([tm(1), tm(end)])
            ylim([10*floor(min(hr_consensus)/10), 10*ceil(max(hr_consensus)/10)])

            if SAVE_FIGURES
                outfname = fullfile(output_results_folder, strcat(full_fname, '_', ch_names{ch}, '_rpeaks.pdf'));
                set(fig, 'Units', 'Inches');
                screenPosition = get(fig, 'Position');
                fig.PaperUnits = 'Inches';
                fig.PaperSize = [screenPosition(3) screenPosition(4)];
                print(fig, outfname, '-dpdf', '-r300')
                % saveas(fig, outfname, '.pdf');

            end
        end

        if SAVE_RESULTS
            outfname = fullfile(output_results_folder, strcat(full_fname, '_', ch_names(ch), '_rpeaks.mat'));
            save(outfname, 'peaks_indexes', 'peak_indexes_consensus');
        end

        % T_ELAPSED_PER_second_in_ms = 1000*toc(TSTART) / (length(data) / fs)

        % SET A BREAK-POINT HERE TO SEE THE FIGURES
        %keyboard % for debugging and visualization only; remove this line during batch runs
        % close(fig)
        disp(strcat('record #', num2str(k), ' fname: ', full_fname, ' channel #', num2str(ch)))
    end
    % end
end