close all; clear; clc
%--------------------------------------------------------------------------
w1 = 0.72; % First stage baseline wander removal window size (in seconds)
w2 = 0.87; % Second stage baseline wander removal window size (in seconds)
BASELINE_REMOVAL_APPROACH = 'BP';% 'BP', 'MDMN'; % baseline wander removal method
fs = 250.0;
%--------------------------------------------------------------------------
datafilepath = '../../../../QTDataBaseMat/';
filelist = dir(fullfile([datafilepath, '**/*.mat'])); 
% RecordName = filelist(k).name;
RecordName = 'sel102m.mat'    %'sel47m.mat'
full_fname = string(RecordName(1:end-4));
datafilename = fullfile(filelist(1).folder, RecordName);
data = load(datafilename);
data = data.val;
size(data)
data = data(:, 1 : round(30*fs)); % select a short segment
size(data)
%------------------
det = 5000;
t = (0 : length(data) - 1)/fs;
figure(1);
subplot(2,1,1);
plot(t(1:det), data(1,1:det), 'linewidth', 1);
legend('Raw ECG ch 1');
subplot(2,1,2); 
plot(t(1:det), data(2,1:det), 'linewidth', 1);
legend('Raw ECG  ch 2');
%--------------------------------------------------------------------------
switch(BASELINE_REMOVAL_APPROACH)
    case 'BP'
        data = data - LPFilter(data, 5.0/fs);
        data = LPFilter(data, 80.0/fs);
    case 'MDMN'
        wlen1 = round(w1 * fs);
        wlen2 = round(w2 * fs);
        for jj = 1 : size(data, 1)
            bl1 = BaseLine1(data(jj, :), wlen1, 'md');
            data(jj, :) = data(jj, :) - BaseLine1(bl1, wlen2, 'mn');
        end
    otherwise
        warning('Unknown method. Bypassing baseline wander removal.');
end
%--------------------------------------------------------------------------
GPfilterparams.bins = 300; % number of phase domain bins
GPfilterparams.BEAT_AVG_METHOD = 'MEAN';%'MEDIAN'; % 'MEAN' or 'MEDIAN'
GPfilterparams.NOISE_VAR_EST_METHOD = 'AVGLOWER'; %'MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE'
GPfilterparams.p = 0.5;
GPfilterparams.avg_bins = 10;
GPfilterparams.wlen_phase = 3;
GPfilterparams.wlen_time = 3;
GPfilterparams.SMOOTH_PHASE = 'BYPASS'; %'MA' or 'GAUSSIAN';
GPfilterparams.gaussianstd = 1.0;
GPfilterparams.plotresults = 0;

GPfilterparams.nvar_factor = 1.0; % noise variance over/under estimation factor (1 by default)
%--------------------------------------------------------------------------
ch = 1;
sig = data(ch, :)/100;
SNR_pre_set = 20;
%--------------------------------------------------------------------------
peak_detector_params.saturate = 1;
peak_detector_params.hist_search_th = 0.9;
peak_detector_params.rpeak_search_wlen = 0.3; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
peak_detector_params.filter_type = 'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER'
[peaks, peaks_indexes] = PeakDetectionProbabilistic(sig, fs, peak_detector_params);

f0 = 1.0; % approximate heart rate (in Hz) used for R-peak detection
[peaks, peaks_indexes] = PeakDetection(sig, f0/fs);                  % peak detection
%--------------------------------------------------------------------------
sd = sqrt(var(sig) / 10^(SNR_pre_set/10));
noise = sd * randn(size(sig));
x = sig + noise;
s_power = mean(sig.^2);
noise_power = mean(noise.^2);
GPfilterparams.nvar = noise_power; % We are assuming that the noise power can be estimated very accurately in practice
%--------------------------------------------------------------------------
[data_posterior_est_phase_based, data_prior_est_phase_based] = ECGPhaseDomainMAPFilter(x, peaks, GPfilterparams);
[data_posterior_est_time_based, data_prior_est_time_based] = ECGTimeDomainMAPFilter(x, peaks, GPfilterparams);
data_wavden = wden(x, 'rigrsure', 's', 'sln', 4, 'sym5');
%--------------------------------------------------------------------------
SNR_pre = 10 * log10(s_power / noise_power);
SNR_prior_phase_based = 10 * log10(s_power / mean((data_prior_est_phase_based - sig).^2));
SNR_posterior_phase_based = 10 * log10(s_power / mean((data_posterior_est_phase_based - sig).^2));
SNR_prior_time_based = 10 * log10(s_power / mean((data_prior_est_time_based - sig).^2));
SNR_posterior_time_based = 10 * log10(s_power / mean((data_posterior_est_time_based - sig).^2));
SNR_wavden = 10 * log10(s_power / mean((data_wavden - sig).^2));
%--------------------------------------------------------------------------
t = (0 : length(x) - 1)/fs;
%------------------
figure(2);
subplot(2,1,1);
plot(t(1:det), sig(1:det), 'linewidth', 1);
legend('Clean ECG');
subplot(2,1,2); 
plot(t(1:det), x(1:det), 'linewidth', 1);
legend('Noisy ECG');
%------------------
figure(3);
subplot(2,1,1);
plot(t(1:det), data_prior_est_phase_based(1:det), 'linewidth', 1); 
hold on;
plot(t(1:det), sig(1:det), 'linewidth', 1);
legend(['GP prio (delta SNR =',num2str(SNR_prior_phase_based-SNR_pre),')'], ['Clean ECG (power = ',num2str(mean(sig.^2)),')']);
subplot(2,1,2);
plot(t(1:det), data_prior_est_phase_based(1:det)-sig(1:det), 'linewidth', 1); 
legend(['noise prio (power = ',num2str(mean((data_prior_est_phase_based - sig).^2)),')']);
%------------------
figure(4);
subplot(2,1,1);
plot(t(1:det), data_posterior_est_phase_based(1:det), 'linewidth', 1); 
hold on;
plot(t(1:det), sig(1:det), 'linewidth', 1);
legend(['GP post (delta SNR =',num2str(SNR_posterior_phase_based-SNR_pre),')'], ['Clean ECG (power = ',num2str(mean(sig.^2)),')']);
subplot(2,1,2);
plot(t(1:det), data_posterior_est_phase_based(1:det)-sig(1:det), 'linewidth', 1); 
legend(['noise post (power = ',num2str(mean((data_posterior_est_phase_based - sig).^2)),')']);
%------------------
figure(5);
subplot(2,1,1);
plot(t(1:det), data_wavden(1:det), 'linewidth', 1); 
hold on;
plot(t(1:det), sig(1:det), 'linewidth', 1);
legend(['Wave (delta SNR =',num2str(SNR_wavden-SNR_pre),')'], ['Clean ECG (power = ',num2str(mean(sig.^2)),')']);
subplot(2,1,2);
plot(t(1:det), data_wavden(1:det)-sig(1:det), 'linewidth', 1); 
legend(['noise wave (power = ',num2str(mean((data_wavden - sig).^2)),')']);
%--------------------------------------------------------------------------










% 
% 
% 
% % Process
% if 1
%     % Test script for an ECG denoiser based on a data driven MAP estimator in
%     % the phase and time domains and QT-interval estimation
% 
%     close all
%     clear
%     clc
% 
%     % InPath = './../../../DataFiles/physionet.org/files/qtdb/1.0.0';
%     % OutPath = './../Results';
%     QTResultsFile = './../../../ResultsReza/QTAnalysisResultsOptNoiseVar.csv';
%     %       QTResultsFile = './QTAnalysisResultsOptNoiseVarFullCovMats.csv';
%     % Load data
%     % datafilepath = '../../../DataFiles/Physionet.org/files/ptbdb/1.0.0/'; % Change this path to where you have the .mat data files
%     % fs = 1000.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
% 
%     %datafilepath = './../../../DataFiles/Physionet.org/files/qtdb/1.0.0/'; 
%     datafilepath = './../../../QTDataBase/';
%     % Change this path to where you have the .mat data files
%     fs = 250.0; % Sampling frequency of the data (put it in the loop and read it from the data if not fixed across all records)
% 
%     SNR_List = -5 : 5 : 30;
%     filelist = dir(fullfile([datafilepath, '**/*.mat']));  % get list of all mat files of interest
%     % Make N by 2 matrix of fieldname + value type
%     variable_names_types = [["FileID", "string"]; ...
%         ["Channel", "int32"]; ...
%         ["SignalLen", "int32"]; ...
%         ["NumRPeaks", "int32"]; ...
%         ["Run", "int32"]; ...
%         ["SetInSNR", "double"]; ...
%         ["ActInSNR", "double"]; ...
%         ["SignalVar", "double"]; ...
%         ["NoiseVar", "double"]; ...
%         ["OutSNR_PriorPhaseBased", "double"]; ...
%         ["OutSNR_PostPhaseBased", "double"]; ...
%         ["OutSNR_PriorTimeBased", "double"]; ...
%         ["OutSNR_PostTimeBased", "double"]];
% 
%     % Make table using fieldnames & value types from above
%     ResultsTableEmpty = table('Size',[1,size(variable_names_types,1)],...
%         'VariableNames', variable_names_types(:,1),...
%         'VariableTypes', variable_names_types(:,2));
% 
%     ResultsTable = ResultsTableEmpty;
%     % Baseline wander removal filter
%     w1 = 0.72; % First stage baseline wander removal window size (in seconds)
%     w2 = 0.87; % Second stage baseline wander removal window size (in seconds)
%     BASELINE_REMOVAL_APPROACH = 'BP'; %'MDMN'; % baseline wander removal method
%     OneRowSet = false;
%     for k = 1 : length(filelist) % Sweep over all or specific records
%         %if isequal(filelist(k).name, 'sele0612m.mat')
%         RecordName = filelist(k).name;
%         %         full_fname = fullfile(OutPath, string(RecordName(1:end-4)));
%         full_fname = string(RecordName(1:end-4));
%         datafilename = fullfile(filelist(k).folder, RecordName);
%         data = load(datafilename);
%         data = data.val;
%         data = data(:, 1 : round(20*fs)); % select a short segment
% 
%         %     OutFname = fullfile(OutPath, string([RecordName(1:end-4), '_raw.csv']));
%         %     writematrix(data', OutFname);
%         %     OutFname = fullfile(OutPath, string([RecordName(1:end-4), '_raw']));
% 
%         switch(BASELINE_REMOVAL_APPROACH)
%             case 'BP'
%                 data = data - LPFilter(data, 5.0/fs);
%                 data = LPFilter(data, 80.0/fs);
%             case 'MDMN'
%                 wlen1 = round(w1 * fs);
%                 wlen2 = round(w2 * fs);
%                 for jj = 1 : size(data, 1)
%                     bl1 = BaseLine1(data(jj, :), wlen1, 'md');
%                     data(jj, :) = data(jj, :) - BaseLine1(bl1, wlen2, 'mn');
%                 end
%             otherwise
%                 warning('Unknown method. Bypassing baseline wander removal.');
%         end
% 
%         for ch = 1 : size(data, 1) % sweep over all or a single desired channel
% 
%             sig = data(ch, :);
% 
%             %%%sig = sig(ch, 1 : round(3*fs))
% 
%             peak_detector_params.saturate = 1;
%             peak_detector_params.hist_search_th = 0.9;
%             peak_detector_params.rpeak_search_wlen = 0.3; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
%             peak_detector_params.filter_type = 'MULT_MATCHED_FILTER';%'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER'
%             [peaks, peaks_indexes] = PeakDetectionProbabilistic(sig, fs, peak_detector_params);
% 
%             f0 = 1.0; % approximate heart rate (in Hz) used for R-peak detection
%             [peaks, peaks_indexes] = PeakDetection(sig, f0/fs);                  % peak detection
% 
%             %[stacked_events, num_non_zeros] = EventStacker(sig, peaks_indexes, 300, 'normalized');
% 
%             TableRow = ResultsTableEmpty;
%             [TableRow.QT_Raw, TableRow.median_QT_Raw, TableRow.median_rr_Raw] = QT_analysis_matrix_in(sig', fs, 1, 1);
%             for SNR_pre_set = SNR_List % the desired input SNR
%                 for itr = 1 : 5 % number of repeated runs
%                     sd = sqrt(var(sig) / 10^(SNR_pre_set/10));
%                     noise = sd * randn(size(sig));
%                     x = sig + noise;
%                     s_power = mean(sig.^2);
%                     noise_power = mean(noise.^2);
% 
%                     GPfilterparams.bins = 300; % number of phase domain bins
%                     GPfilterparams.BEAT_AVG_METHOD = 'MEAN';%'MEDIAN'; % 'MEAN' or 'MEDIAN'
%                     GPfilterparams.NOISE_VAR_EST_METHOD = 'AVGLOWER'; %'MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE'
%                     GPfilterparams.p = 0.5;
%                     GPfilterparams.avg_bins = 10;
%                     %                     GPfilterparams.wlen_phase = 3;
%                     %                     GPfilterparams.wlen_time = 3;
%                     GPfilterparams.SMOOTH_PHASE = 'GAUSSIAN'; %'MA' or 'GAUSSIAN';
%                     GPfilterparams.gaussianstd = 1.0;
%                     GPfilterparams.plotresults = 0;
% 
%                     GPfilterparams.nvar = noise_power; % We are assuming that the noise power can be estimated very accurately in practice
% 
%                     GPfilterparams.nvar_factor = 1.0; % noise variance over/under estimation factor (1 by default)
%                     [data_posterior_est_phase_based, data_prior_est_phase_based] = ECGPhaseDomainMAPFilter(x, peaks, GPfilterparams);
%                     [data_posterior_est_time_based, data_prior_est_time_based] = ECGTimeDomainMAPFilter(x, peaks, GPfilterparams);
%                     %[data_posterior_est_phase_based_fullcov, data_prior_est_phase_based_fullcov] = ECGPhaseDomainMAPFilterFullCovariances(x, peaks, GPfilterparams);
% 
%                     data_wavden = wden(x, 'rigrsure', 's', 'sln', 4, 'sym5');
% 
%                     SNR_pre = 10 * log10(s_power / noise_power);
%                     SNR_prior_phase_based = 10 * log10(s_power / mean((data_prior_est_phase_based - sig).^2));
%                     SNR_posterior_phase_based = 10 * log10(s_power / mean((data_posterior_est_phase_based - sig).^2));
%                     SNR_prior_time_based = 10 * log10(s_power / mean((data_prior_est_time_based - sig).^2));
%                     SNR_posterior_time_based = 10 * log10(s_power / mean((data_posterior_est_time_based - sig).^2));
%                     %SNR_prior_phase_based_fullcov = 10 * log10(s_power / mean((data_prior_est_phase_based_fullcov - sig).^2));
%                     %SNR_posterior_phase_based_fullcov = 10 * log10(s_power / mean((data_posterior_est_phase_based_fullcov - sig).^2));
%                     SNR_wavden = 10 * log10(s_power / mean((data_wavden - sig).^2));
% 
%                     % Log results
%                     TableRow.FileID = string(RecordName);
%                     TableRow.Channel = ch;
%                     TableRow.SignalLen = length(sig);
%                     TableRow.NumRPeaks = length(peaks_indexes);
%                     TableRow.mean_rr0 = mean(diff(peaks_indexes))/fs;
%                     TableRow.median_rr0 = median(diff(peaks_indexes))/fs;
%                     TableRow.Run = itr;
%                     TableRow.SetInSNR = SNR_pre_set;
%                     TableRow.ActInSNR = SNR_pre;
%                     TableRow.SignalVar = s_power;
%                     TableRow.NoiseVar = noise_power;
%                     TableRow.OutSNR_PriorPhaseBased = SNR_prior_phase_based;
%                     TableRow.OutSNR_PostPhaseBased = SNR_posterior_phase_based;
%                     TableRow.OutSNR_PriorTimeBased = SNR_prior_time_based;
%                     TableRow.OutSNR_PostTimeBased = SNR_posterior_time_based;
%                     %TableRow.OutSNR_PriorPhaseBasedFullCov = SNR_prior_phase_based_fullcov;
%                     %TableRow.OutSNR_PostPhaseBasedFullCov = SNR_posterior_phase_based_fullcov;
%                     TableRow.OutSNR_WaveDen = SNR_wavden;
% 
%                     %TableRow.OutFname = strcat(full_fname, ',prior_phase,', num2str(ch), ',', num2str(SNR_pre_set), ',', num2str(itr));
%                     [TableRow.QT_PriorPhaseBased, TableRow.median_QT_PriorPhaseBased, TableRow.median_rr_PriorPhaseBased] = QT_analysis_matrix_in(data_prior_est_phase_based', fs, 1, 1);
% 
%                     [TableRow.QT_PostPhaseBased, TableRow.median_QT_PostPhaseBased, TableRow.median_rr_PostPhaseBased] = QT_analysis_matrix_in(data_posterior_est_phase_based', fs, 1, 1);
% 
%                     [TableRow.QT_PriorTimeBased, TableRow.median_QT_PriorTimeBased, TableRow.median_rr_PriorTimeBased] = QT_analysis_matrix_in(data_prior_est_time_based', fs, 1, 1);
% 
%                     [TableRow.QT_PostTimeBased, TableRow.median_QT_PostTimeBased, TableRow.median_rr_PostTimeBased] = QT_analysis_matrix_in(data_posterior_est_time_based', fs, 1, 1);
% 
%                     %[TableRow.QT_PriorPhaseBasedFullCov, TableRow.median_QT_PriorPhaseBasedFullCov, TableRow.median_rr_PriorPhaseBasedFullCov] = QT_analysis_matrix_in(data_prior_est_phase_based_fullcov', fs, 1, 1);
% 
%                     %[TableRow.QT_PostPhaseBasedFullCov, TableRow.median_QT_PostPhaseBasedFullCov, TableRow.median_rr_PostPhaseBasedFullCov] = QT_analysis_matrix_in(data_posterior_est_phase_based_fullcov', fs, 1, 1);
% 
%                     [TableRow.QT_WaveDen, TableRow.median_QT_WaveDen, TableRow.median_rr_WaveDen] = QT_analysis_matrix_in(data_wavden', fs, 1, 1);
% 
%                     % % ml framework
%                     % soi.q = [-0.045; -0.015];
%                     % soi.t = [.1; .5];
%                     % [mlGaussParams, rPeaks, soi,  waveParams, qtInt] = qtParamsGausFit(x', fs, peaks_indexes, [], soi);
%                     %
%                     % % Bys framework
%                     % PrMu.q = mean(mlGaussParams.q(:, :, 1),2);
%                     % PrCov.q = cov(mlGaussParams.q(:, :, 1)');
%                     %
%                     % PrMu.t = mean(mlGaussParams.t(:, :, 1), 2);
%                     % PrCov.t = cov(mlGaussParams.t(:, :, 1)');
%                     % [bysGaussParams, rPeaks, soi, waveParams, qtInt] = qtParamsGausFit(x', fs, peaks_indexes, [], soi, [], [], [], PrMu, PrCov, noise_power);
% 
%                     % Append results
%                     if OneRowSet == false
%                         ResultsTable = TableRow;
%                         OneRowSet = true;
%                     else
%                         ResultsTable = cat(1, ResultsTable, TableRow);
%                     end
% 
%                     if 0
%                         t = (0 : length(x) - 1)/fs;
%                         figure;
%                         plot(t, x);
%                         hold on;
%                         plot(t(peaks_indexes), sig(peaks_indexes),'ro');
%                         grid
%                         xlabel('time (s)');
%                         ylabel('Amplitude');
%                         title('Noisy ECG and the detected R-peaks');
%                         set(gca, 'fontsize', 16)
%                     end
% 
%                     if 0
%                         lg = {};
%                         t = (0 : length(x) - 1)/fs;
%                         figure
%                         plot(t, x); lg = cat(2, lg, 'Noisy ECG');
%                         hold on
%                         plot(t(peaks_indexes), sig(peaks_indexes),'ro', 'markersize', 12); lg = cat(2, lg, 'R-peaks');
%                         plot(t, data_prior_est_time_based, 'linewidth', 2); lg = cat(2, lg, 'ECG prior estimate (time based)');
%                         plot(t, data_posterior_est_time_based, 'linewidth', 2); lg = cat(2, lg, 'ECG posterior estimate (time based)');
%                         plot(t, data_prior_est_phase_based, 'linewidth', 2); lg = cat(2, lg, 'ECG prior estimate (phase based)');
%                         plot(t, data_posterior_est_phase_based, 'linewidth', 2); lg = cat(2, lg, 'ECG posterior estimate (phase based)');
%                         %plot(t, data_prior_est_phase_based_fullcov, 'linewidth', 2); lg = cat(2, lg, 'ECG prior estimate (phase based full-cov.)');
%                         %plot(t, data_posterior_est_phase_based_fullcov, 'linewidth', 2); lg = cat(2, lg, 'ECG posterior estimate (phase based full-cov.)');
%                         plot(t, data_wavden, 'linewidth', 2); lg = cat(2, lg, 'ECG estimate (waveden)');
%                         plot(t, sig, 'linewidth', 2); lg = cat(2, lg, 'Original ECG');
%                         grid
%                         legend(lg)
%                         xlabel('time (s)');
%                         ylabel('Amplitude');
%                         title('Filtering results');
%                         set(gca, 'fontsize', 16);
%                     end
% 
%                     gap = 6.6;
%                     scale = 30;
%                     if 0
%                         t = (0 : length(x) - 1)/fs;
%                         figure
%                         lg = {};
%                         plot(t, sig/scale, 'linewidth', 2); lg = cat(2, lg, 'Original ECG');
%                         hold on
%                         plot(t(peaks_indexes), sig(peaks_indexes)/scale,'ro', 'markersize', 12); lg = cat(2, lg, 'R-peaks');
%                         plot(t, x/scale - gap * 1); lg = cat(2, lg, 'Noisy ECG');
%                         %                             plot(t, sig/scale, 'linewidth', 2); lg = cat(2, lg, 'Original ECG');
%                         %                             plot(t, data_prior_est_time_based, 'linewidth', 2); lg = cat(2, lg, 'ECG prior estimate (time based)');
%                         %                             plot(t, data_posterior_est_time_based, 'linewidth', 2); lg = cat(2, lg, 'ECG posterior estimate (time based)');
%                         plot(t, data_wavden/scale - gap*2, 'linewidth', 2); lg = cat(2, lg, 'Wavelet denoiser');
%                         plot(t, data_prior_est_phase_based/scale - gap*3, 'linewidth', 2); lg = cat(2, lg, 'Proposed (Prior-based)');
%                         plot(t, data_posterior_est_phase_based/scale - gap*4, 'linewidth', 2); lg = cat(2, lg, 'Proposed (Posterior-based)');
%                         grid
%                         %legend(lg)
%                         xlabel('time [s]');
%                         ylabel('Amplitude [mV]');
%                         %title('Filtering results');
%                         set(gca, 'fontsize', 16);
%                         ylim([-30 10])
% 
%                         %                             figure
%                         %                             lg = {};
%                         %                             grid
%                         %                             legend(lg)
%                         %                             xlabel('time (s)');
%                         %                             ylabel('Amplitude(scaled)');
%                         %                             title('Filtering results');
%                         %                             set(gca, 'fontsize', 16);
%                     end
% 
%                     disp('Filtering performance:');
%                     disp([' Input SNR (desired) = ' num2str(SNR_pre_set), 'dB']);
%                     disp([' Input SNR (actual) = ' num2str(SNR_pre), 'dB']);
%                     
%                     disp([' Output SNR (prior phase-based) = ' num2str(SNR_prior_phase_based), 'dB']);
%                     disp([' Delta SNR (prior phase-based) = ' num2str(SNR_prior_phase_based-SNR_pre), 'dB']);
%                     
%                     %disp([' Output SNR (prior phase-based full-cov) = ' num2str(SNR_prior_phase_based_fullcov), 'dB']);
%                     disp([' Output SNR (prior time-based) = ' num2str(SNR_prior_time_based), 'dB']);
%                     disp([' Delta SNR (prior time-based) = ' num2str(SNR_prior_time_based-SNR_pre), 'dB']);
% 
%                     disp([' Output SNR (posterior phase-based) = ' num2str(SNR_posterior_phase_based), 'dB']);
%                     disp([' Delta SNR (posterior phase-based) = ' num2str(SNR_posterior_phase_based-SNR_pre), 'dB']);
%                     %disp([' Output SNR (posterior phase-based full-cov) = ' num2str(SNR_posterior_phase_based_fullcov), 'dB']);
%                     
%                     disp([' Output SNR (posterior time-based) = ' num2str(SNR_posterior_time_based), 'dB']);
%                     disp([' Delta SNR (posterior time-based) = ' num2str(SNR_posterior_time_based-SNR_pre), 'dB']);
% 
%                     disp([' Output SNR (waveden) = ' num2str(SNR_wavden), 'dB']);
%                     disp([' Delta SNR (waveden) = ' num2str(SNR_wavden-SNR_pre), 'dB']);
%                 end
%             end
%         end
%         %end
%     end
%     writetable(ResultsTable, QTResultsFile);
% end
% 
% % Analyze
% % if 1
% % 
% %     clear;
% %     close all;
% %     clc
% % 
% %         QTResultsFile = './QTAnalysisResultsOptNoiseVar.csv';
% %     QTResultsFile = './QTAnalysisResultsOptNoiseVarFullCovMats.csv';
% %     ResultsTable = readtable(QTResultsFile);
% %     SNR_List = -5 : 5 : 30;
% % 
% %     snr_means = zeros(7, length(SNR_List));
% %     snr_stds = zeros(7, length(SNR_List));
% % 
% %     for k = 1: length(SNR_List)
% %         I = find(ResultsTable.SetInSNR == SNR_List(k));
% %         OutSNR_PriorPhaseBased = ResultsTable.OutSNR_PriorPhaseBased(I) - ResultsTable.ActInSNR(I);
% %         OutSNR_PostPhaseBased = ResultsTable.OutSNR_PostPhaseBased(I) - ResultsTable.ActInSNR(I);
% %         OutSNR_PriorTimeBased = ResultsTable.OutSNR_PriorTimeBased(I) - ResultsTable.ActInSNR(I);
% %         OutSNR_PostTimeBased = ResultsTable.OutSNR_PostTimeBased(I) - ResultsTable.ActInSNR(I);
% %         OutSNR_PriorPhaseBasedFullCov = ResultsTable.OutSNR_PriorPhaseBasedFullCov(I) - ResultsTable.ActInSNR(I);
% %         OutSNR_PostPhaseBasedFullCov = ResultsTable.OutSNR_PostPhaseBasedFullCov(I) - ResultsTable.ActInSNR(I);
% %         OutSNR_WaveDen = ResultsTable.OutSNR_WaveDen(I) - ResultsTable.ActInSNR(I);
% % 
% %         snr_means(1, k) = mean(OutSNR_PriorPhaseBased);
% %         snr_means(2, k) = mean(OutSNR_PostPhaseBased);
% %         snr_means(3, k) = mean(OutSNR_PriorTimeBased);
% %         snr_means(4, k) = mean(OutSNR_PostTimeBased);
% %         snr_means(5, k) = mean(OutSNR_PriorPhaseBasedFullCov);
% %         snr_means(6, k) = mean(OutSNR_PostPhaseBasedFullCov);
% %         snr_means(7, k) = mean(OutSNR_WaveDen);
% % 
% %         snr_stds(1, k) = std(OutSNR_PriorPhaseBased);
% %         snr_stds(2, k) = std(OutSNR_PostPhaseBased);
% %         snr_stds(3, k) = std(OutSNR_PriorTimeBased);
% %         snr_stds(4, k) = std(OutSNR_PostTimeBased);
% %         snr_means(5, k) = std(OutSNR_PriorPhaseBasedFullCov);
% %         snr_means(6, k) = std(OutSNR_PostPhaseBasedFullCov);
% %         snr_stds(7, k) = std(OutSNR_WaveDen);
% % 
% %     end
% % 
% %     figure
% %     subplot(211)
% %     plot(SNR_List, snr_means');
% %     legend({'PriorPhaseBased', 'PostPhaseBased', 'PriorTimeBased', 'PostTimeBased', 'PriorPhaseBasedFullCov', 'PostPhaseBasedFullCov', 'WaveDen'}, 'interpreter', 'none');
% %     xlabel('Input SNR (dB)');
% %     ylabel('SNR Improvement Mean (dB)');
% %     set(gca, 'fontsize', 14);
% %     grid
% % 
% %     subplot(212)
% %     plot(SNR_List, snr_stds');
% %     legend({'PriorPhaseBased', 'PostPhaseBased', 'PriorTimeBased', 'PostTimeBased', 'PriorPhaseBasedFullCov', 'PostPhaseBasedFullCov', 'WaveDen'}, 'interpreter', 'none');
% %     xlabel('Input SNR (dB)');
% %     ylabel('SNR Improvement STD (dB)');
% %     set(gca, 'fontsize', 14);
% %     grid
% % 
% % end