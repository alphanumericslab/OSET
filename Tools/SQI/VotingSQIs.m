function [sortCH, sortSQI] = VotingSQIs(SQIoutputs)
% Voting between different SQIs
%
% Fahimeh Jamshidian Tehrani
% February 2020

CH = size(SQIoutputs, 2);

for ch = 1 : CH
        [r c] = find(SQIoutputs' == ch);
        SQI(ch) = sum(r)/6;
end

[sortSQI, sortCH] = sort(SQI);
