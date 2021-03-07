function fetalHR = FetalHeartRateParams(fetalRPeakTimes, fHRbaselineWinLen, fetalBradycardiaTh, fetalTachycardiaTh)
%
% FetalHeartRateParams(fetalRPeakTimes, Fs, fHRbaselineWinLen)
%
% An implementation of the fetal heart rate parameters as defined in:
% Macones, George A., et al. "The 2008 National Institute of Child Health
%   and Human Development workshop report on electronic fetal monitoring:
%   update on definitions, interpretation, and research guidelines." Journal
%   of Obstetric, Gynecologic, & Neonatal Nursing 37.5 (2008): 510-515.
% DOI: https://doi.org/10.1111/j.1552-6909.2008.00284.x
%
% Inputs:
% fetalRPeakTimes: A time series of the fetal R-peak times instants in seconds
% fHRbaselineWinLen: Fetal baseline HR window length in seconds
% fetalBradycardiaTh: fetal bradycardia threshold in BPM
% fetalTachycardiaTh: fetal tachycardia threshold in BPM
%
% Outputs:
% fetalHR: fetal HR in beats per minute
%
% The Open Source Electrophysiological Toolbox, version 3.14, Oct 2020
% Copyright (C) 2020  Reza Sameni
% reza.sameni@gmail.com
%

fetalHRlen = length(fetalRPeakTimes) - 1;
RRIntervals = diff(fetalRPeakTimes);
fetalHR = 60.0./RRIntervals;

fetalHRbaseline = zeros(1, fetalHRlen);
fetalHRVtype = zeros(1, fetalHRlen);

stop = find(fetalRPeakTimes >= fHRbaselineWinLen, 'first');
if(isempty(stop))
    stop = fetalHRlen;
end
fetalHRbaseline(1) = median(fetalHR(1 : stop));
for itr = 1 : 3
    for k = 1 : fetalHRlen
        start = find(fetalRPeakTimes(k + 1) - fHRbaselineWinLen >= fetalRPeakTimes, 'first');
        if(isempty(start))
            start = 1;
        end
        
        stop = find(fetalRPeakTimes(k + 1) + fHRbaselineWinLen <= fetalRPeakTimes, 'last');
        if(isempty(stop))
            stop = fetalHRlen;
        end
        
        fetalHRsegment = fetalHR(start : stop);
        fetalHRbaseline(k) = median(fetalHRsegment);
    end
    
    for k = 1 : fetalHRlen
        var = abs(fetalHR(k) - fetalHRbaseline(k));
        if(var <= 1.0)
            fetalHRVtype(k) = 0; % fetal HRV absent
        elseif(var > 1.0 && var <= 5.0)
            fetalHRVtype(k) = 1; % fetal HRV minimal
        elseif(var > 5.0 && var <= 25.0)
            fetalHRVtype(k) = 2; % fetal HRV moderate
        elseif(var > 25.0)
            fetalHRVtype(k) = 3; % fetal HRV marked
        end
    end
end


% % % for itr = 1 : 3
% % %     for k = 1 : fetalHRlen
% % %         start = find(fetalRPeakTimes(k + 1) - fHRbaselineWinLen >= fetalRPeakTimes, 'first');
% % %         if(isempty(start))
% % %             start = 1;
% % %         end
% % %
% % %         stop = find(fetalRPeakTimes(k + 1) + fHRbaselineWinLen <= fetalRPeakTimes, 'last');
% % %         if(isempty(stop))
% % %             stop = fetalHRlen;
% % %         end
% % %
% % %         fetalHRsegment = fetalHR(start : stop);
% % %         fetalHRbaseline(k) = median(fetalHRsegment);
% % %     end
% % %
% % %     for k = 1 : fetalHRlen
% % %         var = abs(fetalHR(k) - fetalHRbaseline(k));
% % %         if(var <= 1.0)
% % %             fetalHRVtype(k) = 0; % fetal HRV absent
% % %         elseif(var > 1.0 && var <= 5.0)
% % %             fetalHRVtype(k) = 1; % fetal HRV minimal
% % %         elseif(var > 5.0 && var <= 25.0)
% % %             fetalHRVtype(k) = 2; % fetal HRV moderate
% % %         elseif(var > 25.0)
% % %             fetalHRVtype(k) = 3; % fetal HRV marked
% % %         end
% % %     end
% % % end