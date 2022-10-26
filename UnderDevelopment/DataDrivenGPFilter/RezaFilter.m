% Matlab function that for some levels of noise and a number of repetition
% is reading the csv noisy signals (used for GP-filitering) and filters them
% via Reza's matlab GP implementation and produces the corresponding csv files
% containing the filtered signals corresponding to Reza's
% method/implementation.
clear all; close all; clc
%--------------------------------------------------------------------------
% The sampling frequency 
fs = 250;
%--------------------------------------------------------------------------
% The peak detector's parameters as per Reza's implementation
peak_detector_params.saturate = 1;
peak_detector_params.hist_search_th = 0.9;
peak_detector_params.rpeak_search_wlen = 0.3; % MAX detectable HR (in BPM) = 60/rpeak_search_wlen
peak_detector_params.filter_type = 'BANDPASS_FILTER'; %'BANDPASS_FILTER', 'MATCHED_FILTER', 'MULT_MATCHED_FILTER'
%--------------------------------------------------------------------------
% The GP filter's parameters as per Reza's implementation
GPfilterparams.bins = 300; % number of phase domain bins
GPfilterparams.BEAT_AVG_METHOD = 'MEAN';%'MEDIAN'; % 'MEAN' or 'MEDIAN'
GPfilterparams.NOISE_VAR_EST_METHOD = 'AVGLOWER'; %'MIN', 'AVGLOWER', 'MEDLOWER', 'PERCENTILE'
GPfilterparams.p = 0.5;
GPfilterparams.avg_bins = 10;
%                     GPfilterparams.wlen_phase = 3;
%                     GPfilterparams.wlen_time = 3;
GPfilterparams.SMOOTH_PHASE = 'GAUSSIAN'; %'MA' or 'GAUSSIAN' or 'BYPASS';
GPfilterparams.gaussianstd = 1.0;
GPfilterparams.wlen_phase = 3;
GPfilterparams.wlen_time = 3;
GPfilterparams.plotresults = 0;
%------------------------------------
% The noise variance parameters as per Reza's implementation
GPfilterparams.nvar_factor = 1.0; % noise variance over/under estimation factor (1 by default)
%------------------------------------
noise_levels = [30, 25, 20, 15, 10, 5, 0, -5];
reps = [1, 2, 3, 4, 5];
for noise_level = noise_levels
    for rep = reps
        name_csv_input_clean = ['../../../ResultsComparison120sLP/CleanSignals' ,num2str(rep), '.csv'];
        name_csv_input_noisy = ['../../../ResultsComparison120sLP/NoisySignals', num2str(noise_level), 'db_' ,num2str(rep), '.csv'];        
        name_csv_output_post = ['../../../ResultsComparison120sLP/RezaFiltPost', GPfilterparams.SMOOTH_PHASE, num2str(noise_level),'db_' ,num2str(rep), '.csv'];
        name_csv_output_prio = ['../../../ResultsComparison120sLP/RezaFiltPrio', GPfilterparams.SMOOTH_PHASE, num2str(noise_level),'db_' ,num2str(rep), '.csv'];        
        %------------------------------------
        fprintf(['reading clean signals csv from ', name_csv_input_clean, '...\n'])
        CleanSignals = readtable(name_csv_input_clean, 'NumHeaderLines', 0);
        fprintf(['read noisy signals csv from ', name_csv_input_clean, '.\n'])
        %------------------------------------
        fprintf(['reading noisy signals csv from ', name_csv_input_noisy, '...\n'])
        NoisySignals = readtable(name_csv_input_noisy, 'NumHeaderLines', 0);
        [nlin, ncol] = size(NoisySignals);
        fprintf(['read noisy signals csv from ', name_csv_input_noisy, '.\n'])
        %------------------------------------
        % Reza implementation filtering signals one by one
        fprintf(['Reza implementation filtering noisy signals csv for ',num2str(noise_level), 'dB, rep ',num2str(rep), '...\n'])
        %------------------------------------
        REZAFiltPrio = NaN(nlin,ncol);
        REZAFiltPrio(1:nlin,1) = 0:nlin-1;
        %------------------------------------        
        REZAFiltPost = NaN(nlin,ncol);        
        REZAFiltPost(1:nlin,1) = 0:nlin-1;
        %------------------------------------        
        for i = 2:ncol
            %fprintf(['filtering sig' num2str(i-1), '.\n'])
            CleanSignal = table2array(CleanSignals(:,i));
            CleanSignal = CleanSignal(~isnan(CleanSignal));
            %------------------------------------
            NoisySignal = table2array(NoisySignals(:,i));
            NoisySignal = NoisySignal(~isnan(NoisySignal));
            %------------------------------------
            NoiseSignal = NoisySignal - CleanSignal;
            noise_power = mean(NoiseSignal.^2);
            %------------------------------------
            % the nvar is set as the REAL one.
            GPfilterparams.nvar = noise_power; % We are assuming that the noise power can be estimated very accurately
            %------------------------------------
            % Params as per Reza's implementation
            f0 = 1.0;                                       % approximate heart rate (in Hz) used for R-peak detection
            [peaks, peaks_indexes] = PeakDetection(CleanSignal, f0/fs);                  % peak detection
            %------------------------------------
            [RezaFiltPost, RezaFiltPrio, n_var] = ECGPhaseDomainMAPFilter(NoisySignal', peaks, GPfilterparams);
            %[noise_power,n_var] % uncomment to check the fact that the noise variance is the real one 
            REZAFiltPost(1:length(RezaFiltPost),i)= RezaFiltPost;  
            REZAFiltPrio(1:length(RezaFiltPrio),i)= RezaFiltPrio;              
        end
        fprintf(['Reza implementation filtered noisy signals csv for ',num2str(noise_level), 'dB, rep ',num2str(rep), '.\n'])
        %------------------------
        % Export the Reza implementation filtered posterior signals in a csv file
        fprintf(['exporting the the Reza implementation filtered posterior in ', name_csv_output_post, '...\n'])
        writematrix(REZAFiltPost,name_csv_output_post) 
        fprintf(['exported the the Reza implementation filtered posterior in ', name_csv_output_post, '. DONE \n\n'])
        %------------------------
        % Export the Reza implementation filtered prior signals in a csv file
        fprintf(['exporting the the Reza implementation filtered prior in ', name_csv_output_prio, '...\n'])
        writematrix(REZAFiltPrio,name_csv_output_prio) 
        fprintf(['exported the the Reza implementation filtered prior in ', name_csv_output_prio, '. DONE \n\n'])        
    end
end