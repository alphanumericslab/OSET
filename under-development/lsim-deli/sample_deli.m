

clear;
close all;
clc;


oset_path = 'D:\projects\toolboxes\OSET'; % enter the path of OSET toolbox
addpath(genpath(oset_path))


db_folder = pwd; % if you run from lsim-deli folder otherwise
% db_folder = 'D:\projects\toolboxes\OSET\under-development\lsim-deli'
local_db_files = dir([db_folder '/*.mat']); % list of all mat files


%%

for m = 1:length(local_db_files)

    clc
    close all

    disp(m)
    toc
    tic
    % load data and rpeaks
    in_fname = local_db_files(m).name(1:end-4);

    load([db_folder,'/',in_fname,'.mat']);
    index_R =  true_position.R;

    % ECG preprocessing
    ecg_denoised = ecg(:,1)';

    % NOTCH FILTERING THE ECG
    fc = 50.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_denoised = filtfilt(b, a, ecg_denoised); % zero-phase non-causal filtering
    ecg_denoised = ecg_denoised - movmean(movmedian(ecg_denoised,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);
    ecg_denoised = lp_filter_zero_phase(ecg_denoised, 30/fs);


    % R-peak detection
    peak_detector_params.RETURN_SIGNAL_PEAKS = true; % return signal peaks or energy envelope peaks
    peak_detector_params.PLOT_RESULTS = false; % plot the results using the internal plot function of peak_det_likelihood or not
    peak_detector_params.PLOT_DIAGNOSTIC = false; % diagnostic mode (do not activate unless diving deep into the code! run only on short segments, since many figures are created)
    peak_detector_params.verbose = false; % reports all the default values for the internal parameters of peak_det_likelihood, which can be modified through this data structure if needed.
    peak_detector_params.REFINE_PEAKS = true;
    overlap_time = 1.0; % overlap between segments for continuity (1.0-2.0 seconds is enough)
    seg_len_time = 10.0; % segment length in seconds

    [peaks, ecg_rpeaks_index, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(ecg_denoised, fs, seg_len_time, overlap_time, peak_detector_params);
    ecg_rpeaks_index = peak_indexes_consensus;

    % LSIM-Deli
    flag_post_processing = 1;
    lsim_positions = fiducial_det_lsim( ecg_denoised, ecg_rpeaks_index, fs, flag_post_processing);


    % plot the results
    ecg_rpeaks_index_p = lsim_positions.R;
    ecg_rpeaks_index_p(isnan(ecg_rpeaks_index_p))=[];

    ecg_QRSon_index = lsim_positions.QRSon; ecg_QRSon_index_p = ecg_QRSon_index;  ecg_QRSon_index_p(isnan(ecg_QRSon_index_p))=[];
    ecg_QRSon_true = true_position.QRSon; ecg_QRSon_wavdet_p = ecg_QRSon_true;  ecg_QRSon_wavdet_p(isnan(ecg_QRSon_wavdet_p))=[];

    ecg_QRSoff_index = lsim_positions.QRSoff; ecg_QRSoff_index_p = ecg_QRSoff_index;  ecg_QRSoff_index_p(isnan(ecg_QRSoff_index_p))=[];
    ecg_QRSoff_true = true_position.QRSoff; ecg_QRSoff_wavedet_p = ecg_QRSoff_true;  ecg_QRSoff_wavedet_p(isnan(ecg_QRSoff_wavedet_p))=[];

    ecg_Pon_index = lsim_positions.Pon;  ecg_Pon_index_p = ecg_Pon_index; ecg_Pon_index_p(isnan(ecg_Pon_index_p))=[];
    ecg_P_index = lsim_positions.P;  ecg_P_index_p = ecg_P_index; ecg_P_index_p(isnan(ecg_P_index_p))=[];
    ecg_Poff_index = lsim_positions.Poff;  ecg_Poff_index_p = ecg_Poff_index; ecg_Poff_index_p(isnan(ecg_Poff_index_p))=[];

    ecg_Ton_index = lsim_positions.Ton';     ecg_Ton_index_p = ecg_Ton_index;  ecg_Ton_index_p(isnan(ecg_Ton_index_p))=[];
    ecg_Ton_true = true_position.Ton';     ecg_Ton_wavedet_p = ecg_Ton_true;  ecg_Ton_wavedet_p(isnan(ecg_Ton_wavedet_p))=[];

    ecg_T_index = lsim_positions.T';     ecg_T_index_p = ecg_T_index;  ecg_T_index_p(isnan(ecg_T_index_p))=[];
    ecg_T_true = true_position.T';     ecg_T_wavedet_p = ecg_T_true;  ecg_T_wavedet_p(isnan(ecg_T_wavedet_p))=[];

    ecg_Toff_index = lsim_positions.Toff';  ecg_Toff_index_p = ecg_Toff_index;  ecg_Toff_index_p(isnan(ecg_Toff_index_p))=[];
    ecg_Toff_true = true_position.Toff';  ecg_Toff_wavedet_p = ecg_Toff_true;  ecg_Toff_wavedet_p(isnan(ecg_Toff_wavedet_p))=[];


    a = figure('Position', [130 130 1500 800]);
    lg = {};
    plot(t_second,ecg_denoised,LineWidth=1.5) ;lg = cat(2, lg, {'ECG'});
    hold on
    ecg_plot = ecg_denoised;
    plot(t_second(ecg_rpeaks_index_p),ecg_plot(ecg_rpeaks_index_p),'*',MarkerSize=6,LineWidth=1);lg = cat(2, lg, {'R'});
    plot(t_second(ecg_Ton_index_p),ecg_plot(ecg_Ton_index_p),'m+',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'Ton'});
    plot(t_second(ecg_T_index_p),ecg_plot(ecg_T_index_p),'m*',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'T'});
    plot(t_second(ecg_Toff_index_p),ecg_plot(ecg_Toff_index_p),'mx',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'Toff'});

    if ~isempty(ecg_Ton_wavedet_p)
        plot(t_second(ecg_Ton_wavedet_p),ecg_plot(ecg_Ton_wavedet_p),'g+',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-Ton'});
    end
    if ~isempty(ecg_T_wavedet_p)
        plot(t_second(ecg_T_wavedet_p),ecg_plot(ecg_T_wavedet_p),'g*',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-T'});
    end
    if ~isempty(ecg_Toff_wavedet_p)
        plot(t_second(ecg_Toff_wavedet_p),ecg_plot(ecg_Toff_wavedet_p),'gx',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-Toff'});
    end
    plot(t_second(ecg_Pon_index_p),ecg_plot(ecg_Pon_index_p),'k+',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'Pon'});
    plot(t_second(ecg_P_index_p),ecg_plot(ecg_P_index_p),'k*',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'P'});
    plot(t_second(ecg_Poff_index_p),ecg_plot(ecg_Poff_index_p),'kx',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'Poff'});

    plot(t_second(ecg_QRSon_index_p),ecg_plot(ecg_QRSon_index_p),'rx',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'QRSon'});
    plot(t_second(ecg_QRSon_wavdet_p),ecg_plot(ecg_QRSon_wavdet_p),'go',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-QRSon'});

    plot(t_second(ecg_QRSoff_index_p),ecg_plot(ecg_QRSoff_index_p),'rx',MarkerSize=12,LineWidth=2);lg = cat(2, lg, {'QRSoff'});
    plot(t_second(ecg_QRSoff_wavedet_p),ecg_plot(ecg_QRSoff_wavedet_p),'go',MarkerSize=10,LineWidth=2);lg = cat(2, lg, {'EXPERT-QRSoff'});

    grid on
    legend(lg,'Interpreter' ,'latex','orientation','horizontal','FontSize',14)
    xlabel('time (sec)',Interpreter='latex',FontSize=14)


end



