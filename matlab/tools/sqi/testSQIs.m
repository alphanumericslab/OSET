% A test code for Signal quality indexes based on different methods
% The channel selection method tries to rank the components according to
% their similarity to the fetal mixtures.
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% The Open Source Electrophysiological Toolbox, version 3.14, February 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/





clc
clear all
close all

% Data loading
filename = 'ExtractfECGecgca771SNR10.mat';
load(filename);

% Parameter nitializing
[CH, ~] = size(fECGdata);
fs = 1000;
segLen = 10*fs; % segmentation length
noSegLen = length(fECGdata)/segLen; % # of segments
average_peak_detection_rate = 2.5; %average FECG peak detection rate
flgInput = 1; % Asking user for channel ranking; 1:enable, 0:disable
winLen = 120; % window length for SQI5 algorithm
flgOutput = 1; % Print Results
peaks = cell(CH); % R-peak impulse train for each channel
SQIoutputs = cell(noSegLen); % Calculated SQI(s) for each segment of data
visualIND = cell(noSegLen); % Visual SQI for each segment of data
VotedCH = cell(noSegLen); % Voted SQI for each segment of data

fig = figure(1); set(gcf,'Visible', 'off');
% Applying SQIs for each segment of data
for seg = 1 : noSegLen
    
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
    [valindSQI1, indSQI1] = SQI1(JADEcomp, peaks);
    
    [valindSQI2, indSQI2] = SQI2(JADEcomp, peaks);
    
    [valindSQI3, indSQI3] = SQI3(JADEcomp);
    
    [valindSQI4, indSQI4] = SQI4(JADEcomp);
    
    [valindSQI5, indSQI5] = SQI5(JADEcomp, winLen);
    
    [valindSQI6, indSQI6] = SQI6(JADEcomp, peaks);
    
    SQIoutputs{seg} = [indSQI1; indSQI2; indSQI3; indSQI4; indSQI5; indSQI6;];

    close(fig);
    fig = PlotECG(JADEcomp,CH,'b',fs);
    set(gcf,'Visible', 'on');
    if flgInput  % Asking user for channel ranking
        for ch = 1 : CH
            s = input(strcat('Please enetr the channel # of rank',num2str(ch),': '));
            while isempty(s) || s > CH
                s = input(strcat('Please enetr again the channel # of rank',num2str(ch), ': '));
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


