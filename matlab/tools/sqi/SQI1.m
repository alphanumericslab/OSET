function [valindSQI, indSQI] = SQI1(x, peaks)
%
% [valindSQI, indSQI] = SQI1(x, peaks)
% Channel selection based on Pseudo-periodicity
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
%   valindSQI: values of "Pseudo-periodicity" measure
%   corresponding to the channel number of "indSQI" parameter
%   indSQI: ranked channel number according to the "Pseudo-periodicity" measure
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
mk = cell(CH);
RWAmatrix = cell(CH);
ind = zeros(1, CH);

for ch = 1: CH

    Peakloc{ch} = find(peaks{ch});
    i = 1;
    T = mean(diff(Peakloc{ch}));
    winL = round(T/2);
    while Peakloc{ch}(i)<winL
        RWAmatrix{ch}(i,1:2*winL) = [zeros(1,winL-Peakloc{ch}(i)) x(ch,1:Peakloc{ch}(i)+winL)];
        i = i +1;
    end
    for j = i : length(Peakloc{ch})-1
        if Peakloc{ch}(j)+winL <= length(x(ch, :))
            RWAmatrix{ch}(j, 1:2*winL) = x(ch,Peakloc{ch}(j)-winL+1 : Peakloc{ch}(j)+winL);
        else
            RWAmatrix{ch}(j, 1:2*winL) = [x(ch, Peakloc{ch}(j)-winL+1:end) zeros(1,2*winL-length(x(ch, Peakloc{ch}(j)-winL+1:end))) ];
        end
    end
    if 2*winL > length(x(ch, Peakloc{ch}(end)-winL+1:end))
        RWAmatrix{ch}(length(Peakloc{ch}), 1:2*winL) = [x(ch, Peakloc{ch}(end)-winL+1:end) zeros(1,2*winL-length(x(ch, Peakloc{ch}(end)-winL+1:end))) ];
    else
        RWAmatrix{ch}(length(Peakloc{ch}), 1:2*winL) = x(ch, Peakloc{ch}(end)-winL+1:Peakloc{ch}(end)+winL);
    end
    mk{ch} = RWAverage(RWAmatrix{ch});

    ind(ch) = mean(mean((RWAmatrix{ch} - ones(size(RWAmatrix{ch},1),1)* mk{ch}).^2));
end

% ranking channels according to the "Pseudo-periodicity" measure
[valindSQI, indSQI] = sort(ind);