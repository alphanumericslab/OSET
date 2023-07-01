function [valindSQI, indSQI] = SQI6(x, peaks)
%
% [valindSQI, indSQI] = SQI6(x, peaks)
% Channel selection based on measure of PICA periodicity
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% input:
%   x: input fetal ECG record
%   peaks: R-peak impulse train
% output:
%   valindSQI: values of "PICA periodicity" measure
%   corresponding to the channel number of "indSQI" parameter
%   indSQI: ranked channel number according to the "PICA periodicity" measure
%
% The Open Source Electrophysiological Toolbox, version 3.14, February 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/



% Comment out for individual using
% % %     Remove the mean
% % mn = mean(x,2)*ones(1,size(x,2));
% % x = x - mn;
% % 
% % %     Normalize the variance of x
% % x = x./var(x,0,2);

% Parameter nitializing
[CH, ~] = size(x);
ind = zeros(1, CH);
T0 = cell(CH);
T1 = cell(CH);

for ch = 1: CH
    Peakloc = find(peaks{ch});
    T = mean(diff(Peakloc));
    [T0{ch},T1{ch}] = SynchPhaseTimes2(peaks{ch});
    be4peak = 1: T0{ch}(1)-1;
    T0{ch} = [be4peak T0{ch}];
    T1{ch} = [be4peak+round(T) T1{ch}];
    
    term1 = mean(x(ch, T0{ch}).*x(ch, T1{ch}));
    term2 = mean(x(ch,:).^2);
    ind(ch) = term1 / term2;
end

% ranking channels according to the "PICA periodicity" measure
[valindSQI, indSQI] = sort(ind, 'descend');