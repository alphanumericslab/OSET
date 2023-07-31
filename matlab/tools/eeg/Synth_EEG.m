function eeg_synth = Synth_EEG(sig, fs, varargin)
% 
% eeg_synth = Synth_EEG(sig, fs)
% eeg_synth = Synth_EEG(sig, fs, duration, win, AR_ord)
% 
% *************************************************************************
% * Generating synthetic EEG signal using Autoregressive (AR) model and   *
% * innovation filter. Could be used in cross-validations with real EEG   *
% * signal. Refer to the User Guide for further details.                  *
% *************************************************************************
% 
%             Inputs:
%                       sig : input raw EEG signal
%                       fs : sampling frequency (Hz)
%            (optional) duration : total required signal duration (sec)
%            (optional) win : temporal window length to have a stationary 
%                             signal (Sec)
%            (optional) AR_ord: order of AR model
% 
%             Defaults:
%                       duration = length of input raw EEG signal
%                       win = 3 or 4 or 5 sec (depends on input EEG signal)
%                       AR_ord = 20
% 
%             *NOTE: While specifying a value to one of the parameters
%                    having default values, an empty bracket [] must be
%                    used for non-specified parameters. If you're not using
%                    these parameters, empty bracket is not required.
%             Outputs: 
%                       eeg_synth : a synthetic EEG signal with spectral
%                                   characteristics similar to a real EEG
%                                   signal.
% 
% This program is provided by ESMAEIL SERAJ (esmaeil.seraj09@gmail.com).
% Please make sure to reference BOTH the original studies [1-2] and the 
% OSET [3] to help others find these items.
% 
%     [1] Esmaeil Seraj, Reza Sameni. ”Robust Electroencephalogram Phase 
%         Estimation with Applications in Brain-computer Interface Systems” 
%         Physiological Measurements (2017)
%     [2] Reza Sameni and Esmaeil Seraj, “A Robust Statistical Framework 
%         for Instantaneous Electroencephalogram Phase and Frequency 
%         Analysis” Physiological Measurements (2017)                     
%     [3] R. Sameni, The Open-Source Electrophysiological Toolbox (OSET), 
%         version 3.1 (2014). URL http://www.oset.ir
%         Released under the GNU General Public License
%         Copyright (C) 2012  Reza Sameni
%         Shiraz University, Shiraz, Iran
%         reza.sameni@gmail.com 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.
% 

%%-Checking inputs and assigning default values-%%
sig = sig(:)';
if nargin<2
    error('***Not Enough Input Arguments: Please check the Instructions***')
else
    if nargin~=5
        error('***An Empty Bracket [] Must be Used for the Values NOT Being Specified***')
    end
end

[m1, n1] = size(sig);
if m1~=1
    error('***The input EEG signal must be a SINGLE-channel vector***')
end
if nargin==5
    if isempty(varargin{1})
        duration = n1/fs;
    else
        duration = varargin{1};
    end
    if isempty(varargin{2})
        for i=3:9
            if rem(duration, i)==0
                win = i;
            end
        end
        if isempty(win)
            error('***Please change the total required duration to another close value***')
        end
    else
        win = varargin{2};
    end
    if isempty(varargin{3})
        AR_ord = 20;
    else
        AR_ord = varargin{3};
    end
end
if rem(duration, win)~=0
    error('***The window length has to be dividable by the total required signal duration***')
end
if ~(isinteger(fs))&&(isinteger(duration))
   error('***The values assigned to "fs" and "duration" must be positive integers***') 
end

%---matching the real signal into required duration for synthetic eeg
%signal---%
sig_length = fs*duration;
sig = sig(1, 1:sig_length);
[~ ,n] = size(sig);
win_length = win*fs;

%---generating synthetic eeg using inovation filter---%
AR_coeffs = zeros(duration/win, AR_ord+1);
for i=1:win_length:n
    AR_coeffs(((i-1)/win_length)+1, :) = arburg(sig(1, i:i+(win_length-1)), AR_ord);
end
eeg_synth = zeros(duration/win, win_length);
for j=1:duration/win
    white_n = randn(1, win_length);          % white noise
    eeg_synth(j, :) = filter(1, AR_coeffs(j, :), white_n);
end
eeg_synth = reshape(eeg_synth', 1, []);

end