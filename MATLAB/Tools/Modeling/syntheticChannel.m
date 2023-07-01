function [synCh, beat] = syntheticChannel(x, winL, Peakloc)
%
% [synCh, beat] = syntheticChannel(x, winL, Peakloc)
% synthetic channel construction based on Robust weighted averaging
%
% Fahimeh Jamshidian Tehrani,
% jamshidian.t@gmail.com
% September 2022
%
% input:
%   x: input ECG record it can be multichannel
%   winL: Half of the R-R interval
%   Peakloc: peak indexes
% output:
%   synCh: Synthetic channel
%   beat: Synthetic beat
%
% The Open Source Electrophysiological Toolbox, 
% Released under the GNU General Public License
% https://gitlab.com/rsameni/OSET/


[CH, N] = size(x);
i = 1; % strat from the first beat

for chRWA = 1: CH
    if Peakloc(1)<winL % if the first beat is not complete 
        RWAmatrix(1,1:2*winL) = [zeros(1,winL-Peakloc(1)) x(chRWA,1:Peakloc(1)+winL)];
        i = 2; % continue for the second beat
    end
    for j = i : length(Peakloc)-1 % for all beats
        RWAmatrix(j,1:2*winL) = x(chRWA,Peakloc(j)-winL+1 : Peakloc(j)+winL);
    end
    beat_tmp(chRWA,:) = RWAverage(RWAmatrix);
    
end
% average over all beats
beat = mean(beat_tmp, 1);

% synthetic channel construction
% repeat the beat at every peak location
synCh = zeros(size(beat,1), N);
if Peakloc(1)<winL % if the first beat is not complete
    synCh(:,1:Peakloc(1)+winL+1) = beat(:,1:Peakloc(1)+winL+1);
    i = 2; % continue for the second beat
end
for j = i : length(Peakloc)-1
    synCh(:,Peakloc(j)-winL+1:Peakloc(j)+winL) = beat;
end
if Peakloc(end)+winL > size(synCh,2) % if the last beat is not complete
    synCh(:,Peakloc(end)-winL+1:end) = beat(:,1:end - (Peakloc(end)+winL- size(synCh,2)));
else
    synCh(:,Peakloc(end)-winL+1:Peakloc(end)+winL) = beat;
end

end