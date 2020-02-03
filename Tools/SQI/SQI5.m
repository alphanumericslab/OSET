function [valindSQI, indSQI] = SQI5(x, winLen)
%
% [valindSQI, indSQI] = SQI5(x, winLen)
% Channel selection based on channel power in sliding window
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% input:
%   x: input JADE component of fetal ECG record
%   winLen: sliding window length
% output:
%   valindSQI: values of "channel power in sliding window" measure
%   corresponding to the channel number of "indSQI" parameter
%   indSQI: ranked channel number according to the "channel power in sliding window" measure
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
yk = zeros(size(x));

for ch = 1: CH
    yk(ch,:) = movmean(x(ch, :).^2, winLen);
    ind(ch) = var(yk(ch,:));
end

% ranking channels according to the "channel power in sliding window" measure
[valindSQI, indSQI] = sort(ind, 'descend');