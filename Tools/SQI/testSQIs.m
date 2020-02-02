% A test code for Signal quality indexes based on different methods
% The channel selection method tries to rank the components according to
% their similarity to the fetal mixtures.
%
% Fahimeh Jamshidian Tehrani
% February 2020


clc
clear all
close all

% Parameter nitializing
filename = 'ExtractfECGecgca771SNR10.mat';
fs = 1000;
segLen = 10*fs; % segmentation length
average_peak_detection_rate = 2.5; %average FECG peak detection rate
flgInput = 1; % Asking user for channel ranking; 1:enable, 0:disable
winLen = 120; % window length for SQI5 algorithm
flgOutput = 1; % Print Results

% Data loading
load(filename);
[CH, ~] = size(fECGdata);

% Applying SQIs for each segment of data
for seg = 1 : length(fECGdata)/segLen
    
    segIND = (seg-1)*segLen+1 : seg*segLen;
    
    % applying JADE
    CC = jadeR(fECGdata(:,segIND), CH);
    JADEcomp = CC*fECGdata(:,segIND);
    
    % Remove the mean
    mn = mean(JADEcomp,2)*ones(1,size(JADEcomp,2));
    JADEcomp = JADEcomp - mn;
    
    % Normalize the variance of x
    JADEcomp = JADEcomp./var(JADEcomp,0,2);
    
    % Peak detection
    for ch = 1:CH
        peaks{ch} = PeakDetection(JADEcomp(ch,:), average_peak_detection_rate/fs, 0);
    end
    
    % Call the signal quality indexe algorithms
    [valindSQI1, indSQI1] = sort(SQI1(JADEcomp, peaks));
    
    [valindSQI2, indSQI2] = sort(SQI2(JADEcomp, peaks));
    
    [valindSQI3, indSQI3] = sort(SQI3(JADEcomp),'descend');
    
    [valindSQI4, indSQI4] = sort(SQI4(JADEcomp),'descend');
    
    [valindSQI5, indSQI5] = sort(SQI5(JADEcomp, winLen),'descend');
    
    [valindSQI6, indSQI6] = sort(SQI6(JADEcomp, peaks),'descend');
    
    SQIoutputs{seg} = [indSQI1; indSQI2; indSQI3; indSQI4; indSQI5; indSQI6;];
    
    close all
    figure(1);
    PlotECG(JADEcomp,CH,'b',fs);
    if flgInput  % Asking user for channel ranking
        for ch = 1 : CH
            s = input(strcat('Please enetr channel of rank',num2str(ch),': '));
            while isempty(s) || s > 10
                s = input(strcat('Please enetr again the channel of rank',num2str(ch), ': '));
            end
            visualIND{seg}(ch) = s;
        end
    end
    
    % Voting
    [VotedCH{seg}, sortSQI] = VotingSQIs(SQIoutputs{seg});
    
    % Output results
    if flgOutput
        disp('Results of 1 -> 6 SQIs');
        disp(SQIoutputs{seg});
        disp('Result of Voted SQI');
        disp(VotedCH{seg});
        if flgInput
            disp('Result of Visual SQI');
            disp(visualIND{seg});
        end
    end
end


