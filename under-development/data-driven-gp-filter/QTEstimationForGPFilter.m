% Matlab function that is applying the QT detector over signals retrieved
% in CSV form and saves the results in CSV as well.


clear all; close all; clc
%------------------------
%method = 'WaveFilt'
% method = 'GPFiltPrio'
method = 'GPFiltPost'
%------------------------
% Read the GP filtered signals from the csv file
% snr_level = 'Raw';
% snr_level = 30;
% rep = 5;

snr_levels = [30, 25, 20, 15, 10, 5, 0, -5];
% snr_levels = ['Raw'];
reps = [1, 2, 3, 4, 5];
for snr_level = snr_levels
    for rep = reps

        if(ischar(snr_level))
            name_csv_input  = ['../../../../ResultsComparison120s/CleanSignals', num2str(rep)];
            name_csv_output = ['../../../../ResultsComparison120s/QTMedians', snr_level, 'db_', num2str(rep), '.csv'];
        else
            name_csv_input  = ['../../../../ResultsComparison120s/', method, num2str(snr_level), 'db_', num2str(rep), '.csv'];
            name_csv_output = ['../../../../ResultsComparison120s/', 'QTMedians', method, num2str(snr_level), 'db_', num2str(rep), '.csv'];
        end
        %------------------------
        if(ischar(snr_level))
            fprintf('reading the csv file for Raw, REP = %d... \n', rep)
        else
            fprintf('reading the csv file for SNR = %d, REP = %d ... \n', snr_level, rep)
        end
        T = readtable(name_csv_input, 'NumHeaderLines', 0);
        [nlin, ncol] = size(T);
        fprintf('read the csv in table with dimensions %d and %d. \n', nlin, ncol)
        %------------------------
        % Compute the QT and RR medians & save them in a matrix
        fprintf('computing the QT and RR medians... \n')
        for i = 2:ncol
            GPFilteredSig = T(:,i);
            [QT, RR] = QT_analysis_single_lead(GPFilteredSig.Variables,250);
            %------------------------
            medianQT=nanmedian(QT); % for median results
            medianRR=nanmedian(RR); % for median of rr    
            medians = [medianQT;medianRR];
            %------------------------    
            if i == 2
                MEDIANS = medians;
            else 
                MEDIANS = [MEDIANS, medians];
            end    
        end
        [nlin, ncol] = size(MEDIANS);
        fprintf('computed the QT and RR medians and saved in matrix with dimension %d and %d. \n', nlin, ncol)
        %------------------------
        % Export the QT and RR medians & export them in a csv file
        fprintf('saving the the QT and RR medians in csv file... \n')
        writematrix(MEDIANS,name_csv_output) 
        fprintf('saved the the QT and RR medians in csv file. \n')

    end
end