function [UC, UC_onsets, UC_counts, UC_counts_avg, Tachysystolic] = UterineContractions(UA, input_type, Fs, alpha, UC_iden_wlen, UC_avg_wlen, Tachysystolic_th, plotflag, varargin)
%
% [UC, UC_onsets, UC_counts, UC_counts_avg, Tachysystolic] = UterineContractions(UA, input_type, Fs, alpha, UC_iden_wlen, UC_avg_wlen, Tachysystolic_th, plotflag, EHG_fl, EHG_fu, EHG_env_win_param)
%
% An implementation of the uterine contraction calculation as defined in:
% Macones, George A., et al. "The 2008 National Institute of Child Health 
%   and Human Development workshop report on electronic fetal monitoring:
%   update on definitions, interpretation, and research guidelines." Journal
%   of Obstetric, Gynecologic, & Neonatal Nursing 37.5 (2008): 510-515.
% DOI: https://doi.org/10.1111/j.1552-6909.2008.00284.x
%
% Inputs:
% UA: Uterine activity signal in the form of electrohysterogram (EHG), intrauterine pressure (IUP), or tocogram (TOCO)
% input_type: 'EHG' or 'IUP'
% Fs: Input sampling frequency
% alpha: IUP detection threshold parameter (between -1 and +1). Set to -1
%       for th = IUP_min, 0 for th = IUP_median, and +1 for th = IUP_max
%   Note: Negative values scale between IUP_min and IUP_median and positive
%       values scale between IUP_median and IUP_max
% UC_iden_wlen: Uterine contraction identification window length
% UC_avg_wlen: Uterine contraction averaging window length
% Tachysystolic_th: Threshold above which the uterine activity is potantially tachysystolic
% plotflag: plot (1) the results or not (0)
% EHG_fl: EHG lower frequency
% EHG_fu: EHG upper frequency
% EHG_env_win_param: Window parameter used for calculating IUP from EHG
%
% Outputs:
% UC: Uterine contraction waveform
% UC_onsets: Uterine contraction event onset pulses
% UC_counts: Uterine contraction counts over the specified windows
% UC_counts_avg: Uterine contraction counts averaged over the specified windows
% Tachysystolic: The potentially tachy-systolic time instants (the fetal HR should be additionally checked according to the above reference)
%
% The Open Source Electrophysiological Toolbox, version 3.14, Oct 2020
% Copyright (C) 2020  Reza Sameni
% reza.sameni@gmail.com
%

if(isequal(input_type, 'EHG'))
    EHG_fl = varargin{1};
    EHG_fu = varargin{2};
    EHG_env_win_param = varargin{3};
    wlen = round(EHG_env_win_param * Fs); % envelope detection window length in samples
    y = bandpass(UA, [EHG_fl, EHG_fu], Fs); % Bandpassed version of the EHG
    IUP = sqrt(filtfilt(ones(1, wlen), wlen, y.^2)); % IUP signal
elseif(isequal(input_type, 'IUP') || isequal(input_type, 'TOCO'))
    IUP = UA; % IUP/TOCO is given as input, no EHG available
else
    error('Unknown method');
end

UA_len = length(UA); % UA signal length

IUP(IUP < 0) = 0; % pressure may only be positive

IUP_min = min(IUP); % min IUP
IUP_median = median(IUP); % median IUP
IUP_max = max(IUP); % max IUP

if(alpha >= 0) % set threshold somewhere between IUP_median and IUP_max
    th = (1 - alpha) * IUP_median + alpha * IUP_max;
else % set threshold somewhere between IUP_min and IUP_median
    th = (1 + alpha) * IUP_median - alpha * IUP_min;
end

UC = IUP >= th; % Find uterine contraction samples
UC_onsets = diff([0 UC]) == 1; % Find uterine contraction onsets

% Count the number of uterine contractions over a non-causal sliding window of length UC_iden_wlen
UC_counts = zeros(1, UA_len);
UC_iden_half_wlen = round(UC_iden_wlen * Fs / 2);
for k = 1 : UA_len
    if(k <= UC_iden_half_wlen)
        start = 1;
    else
        start = k - UC_iden_half_wlen;
    end
    if(k < UA_len - UC_iden_half_wlen)
        stop = k + UC_iden_half_wlen;
    else
        stop = UA_len;
    end
    UC_counts(k) = sum(UC_onsets(start:stop));
end

% Average the number of uterine contractions over a non-causal sliding window of length UC_avg_wlen
UC_counts_avg = zeros(1, UA_len);
UC_avg_half_wlen = round(UC_avg_wlen * Fs / 2);
for k = 1 : UA_len
    if(k <= UC_avg_half_wlen)
        start = 1;
    else
        start = k - UC_avg_half_wlen;
    end
    if(k < UA_len - UC_avg_half_wlen)
        stop = k + UC_avg_half_wlen;
    else
        stop = UA_len;
    end
    UC_counts_avg(k) = mean(UC_counts(start:stop));
end

% The potentially tachy-systolic time instants
Tachysystolic = UC_counts_avg > Tachysystolic_th;

% Plot the results if demanded
if(plotflag)
    t = (0 : UA_len - 1)/Fs;
    figure
    hold on
    plot(t, UA);
    plot(t, IUP);
    plot(t, th(ones(1, UA_len)));
    grid
    hold off
end