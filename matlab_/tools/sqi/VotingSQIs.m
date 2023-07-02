function [sortCH, sortSQI] = VotingSQIs(SQIInds)
%
% [sortCH, sortSQI] = VotingSQIs(SQIInds)
% Voting between different SQIs
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@cse.shirazu.ac.ir,
% fahimeh.jt@gmail.com
% February 2020
%
% input:
%   SQIInds: different input SQIs
% output:
%   sortCH: ranked channel number according to Voting between different SQIs
%   sortSQI: score of channels corresponding to the channel number of "sortCH" parameter
%
% The Open Source Electrophysiological Toolbox, version 3.14, February 2020
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/

CH = size(SQIInds, 2);
SQI = zeros(1, CH);

for ch = 1 : CH
        [r, ~] = find(SQIInds' == ch);
        SQI(ch) = sum(r)/6;
end

% ranking channels according to their scores
[sortSQI, sortCH] = sort(SQI);
