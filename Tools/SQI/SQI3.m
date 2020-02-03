function [valindSQI, indSQI] = SQI3(x)
%
% [valindSQI, indSQI] = SQI3(x)
% Channel selection based on measure of non-gaussianity
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% input:
%   x: input JADE component of fetal ECG record
% output:
%   valindSQI: values of "spectral energy ratio" measure
%   corresponding to the channel number of "indSQI" parameter
%   indSQI: ranked channel number according to the "spectral energy ratio" measure
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
k1 = 36/(8*sqrt(3)-9);
k2 = 24/(16*sqrt(3)-27);
ind = zeros(1, CH);


for ch = 1: CH
   term1 = k1 * mean(x(ch,:).*exp(-(x(ch,:).^2)/2))^2; 
   term2 = k2 * (mean(exp(-(x(ch,:).^2)/2))-sqrt(2)/2)^2;
   ind(ch) = abs(term1 + term2);
end

% ranking channels according to the "spectral energy ratio" measure
[valindSQI, indSQI] = sort(ind, 'descend');