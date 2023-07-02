% Matlab function that for some levels of noise and a number of repetition
% is reading the csv noisy signals (used for GP-filitering) and filters them
% via wavelet matlab method and produces the corresponding csv files
% containing the waveled-filtered results. 


clear all; close all; clc
%------------------------------------
noise_levels = [30, 25, 20, 15, 10, 5, 0, -5];
reps = [1, 2, 3, 4, 5];
for noise_level = noise_levels
    for rep = reps
        name_csv_input = ['../../../../ResultsComparison120s/NoisySignals', num2str(noise_level), 'db_' ,num2str(rep), '.csv'];
        name_csv_output = ['../../../../ResultsComparison120s/WaveFilt', num2str(noise_level),'db_' ,num2str(rep), '.csv'];
        %------------------------------------
        fprintf(['reading noisy signals csv from ', name_csv_input, '.\n'])
        NoisySignals = readtable(name_csv_input, 'NumHeaderLines', 0);
        [nlin, ncol] = size(NoisySignals);
        fprintf(['read noisy signals csv from ', name_csv_input, '.\n'])
        %------------------------------------
        % Wavelet filtering signals one by one
        fprintf(['Wavelet filtering noisy signals csv for ',num2str(noise_level), 'dB, rep ',num2str(rep), '.\n'])
        WAVEFiltSignal = NaN(nlin,ncol);
        WAVEFiltSignal(1:nlin,1) = 0:nlin-1;
        for i = 2:ncol
            %fprintf(['filtering sig' num2str(i-1), '.\n'])
            NoisySignal = table2array(NoisySignals(:,i));
            NoisySignal = NoisySignal(~isnan(NoisySignal));
            WaveFiltSignal = wden(NoisySignal, 'rigrsure', 's', 'sln', 4, 'sym5');
            WAVEFiltSignal(1:length(WaveFiltSignal),i)= WaveFiltSignal;  
        end
        fprintf(['Wavelet filtered noisy signals csv for ',num2str(noise_level), 'dB, rep ',num2str(rep), '.\n'])
        %------------------------
        % Export the Wavelet filtered signals in a csv file
        fprintf(['exporting the the Wavelet filtered in ', name_csv_output, '.\n'])
        writematrix(WAVEFiltSignal,name_csv_output) 
        fprintf(['exported the the Wavelet filtered in ', name_csv_output, '. DONE \n\n'])
    end
end
