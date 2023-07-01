function [valindSQI, indSQI] = SQI2(x, peaks)
%
% [valindSQI, indSQI] = SQI2(x, peaks)
% Channel selection based on R-peak amplitude constancy
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% input:
%   x: input JADE component of fetal ECG record
%   peaks: R-peak impulse train
% output:
%   valindSQI: values of "R-peak amplitude constancy" measure
%   corresponding to the channel number of "indSQI" parameter
%   indSQI: ranked channel number according to the "R-peak amplitude constancy" measure
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
Peakloc = cell(CH);
ind = zeros(1, CH);

for ch = 1: CH
    Peakloc{ch} = find(peaks{ch});
    ind(ch) = var(x(ch,Peakloc{ch}));
end

% ranking channels according to the "R-peak amplitude constancy" measure
[valindSQI, indSQI] = sort(ind);