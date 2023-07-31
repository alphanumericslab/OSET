function [emg_quant, synch_emg2, ecg_estimate2_bl, time_vec] = emg_quantification(emg_data, fs, emg_onset_sampl, varargin)
% 
% [emg_quant, synch_emg2, ecg_estimate2_bl, time_vec] = emg_quantification(emg_data, fs, emg_onset_sampl)
% [emg_quant, synch_emg2, ecg_estimate2_bl, time_vec] = emg_quantification(emg_data, fs, emg_onset_sampl, duration)
% 
% *************************************************************************
% * EMG Signal Analysis and Quantification                                *
% *************************************************************************
% 
% Usage:    [emg_quant, synch_emg2, ecg_estimate2_bl, time_vec] = emg_quantification(emg_data, fs, emg_onset_sampl)
%           [emg_quant, synch_emg2, ecg_estimate2_bl, time_vec] = emg_quantification(emg_data, fs, emg_onset_sampl, duration)
% inputs:
%           'emg_data': cell array containing emg signals of all trials
%           'fs': sampling frequency (Hz)
%           'emg_onset_sampl': vector containing movement onset samples of
%                              all trials (output of 'bdf2mat_main.m' function)
%     (opt) 'duration': duration of signal required after onset (default: 
%                       duration = 2 (Seconds))
% outputs:
%           'emg_quant': vector containing quantified emg values
%           'synch_emg2': cell array containing synchronized emg trials
%                         based on movement onset
%           'ecg_estimate2_bl': extracted ECG signal from EMG channels
%           'time_vec': time-vector required for plotting quantified EMG
% Note:
%           an empty bracket [] Must be assigned to not-specified values
%
% This program is provided by ESMAEIL SERAJ. Please make sure to cite BOTH 
% the original studies and the User Manual to help others find these items.
% 
% Authors:
% 			Esmaeil Seraj, Karthiga Mahalingam
% Websites:
%			https://github.com/EsiSeraj/ERP_Connectivity_EMG_Analysis
% 			http://oset.ir/category.php?dir=Tools
% 
%    Copyright (C) <2018>  <ESMAEIL SERAJ (eseraj3@gatech.edu)>
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program. If not, see <http://www.gnu.org/licenses/> or
%    write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth
%    Floor, Boston, MA  02110-1301, USA.
% 

%% Checking inputs and assigning default values
if nargin < 3
    error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 3
    duration = 2;
elseif nargin > 3
    if size(varargin, 2) ~= 1
        error('***wong number of variable input arguments. Only <duration> is supposed to be entered***')
    else
        if isempty(varargin)
            duration = 2;
        else
            duration = varargin;
        end
    end
end
if (~(isscalar(fs) && isscalar(duration)))
    error('***sampling frequency and duration have to be scalars***')
end
if ~iscell(emg_data)
    error('***input EMG data needs to be stored in a cell array. Refer to the manual for more details***')
end
if (isscaler(emg_onset_sampl) || iscell(emg_onset_sampl))
    error('***wrong format for EMG onset samples vector. Refer to manual for more details***')
end

%% frequency filtering
ellip_lpfilt_order2 = 2;             % butterworth filter order
lpcutt_off_freq2 = [1 30];           % butterworth filter cut-off frequency
[lp_numerator_coeff2, lp_denominator_coeff2] = ellip(ellip_lpfilt_order2, 0.1,...
    50, lpcutt_off_freq2/(fs/2), 'bandpass');
% freqz(lp_numerator_coeff2, lp_denominator_coeff2)  % filter visualization
ecg_estimate2 = cell(1, length(emg_data));
ecg_estimate2_bl = cell(1, length(emg_data));
for i=1:length(emg_data)
    ecg_estimate2{i} = filtfilt(lp_numerator_coeff2, lp_denominator_coeff2, emg_data{i});
    ecg_estimate2_bl{i} = BaseLine2(ecg_estimate2{i}, 20, 40, 'mn');
end

%% removing ECG
trigger_time = min(emg_onset_sampl);
post_move_time = duration*fs; % duration required after movement onset in seconds
emg_pure2 = cell(1, length(emg_data));
for i=1:length(emg_data)
    emg_pure2{i} = emg_data{i} - ecg_estimate2_bl{i};
end

%% EMG quantification
for i=1:length(emg_pure2)
    sampl_vec = emg_onset_sampl(i)-trigger_time+1:emg_onset_sampl(i)+post_move_time;
    synch_emg2(i, :) = emg_pure2{i}(sampl_vec); % variable size. DO NOT preallocate
end
emg_avg3 = mean(abs(synch_emg2));
time_vec = (1:length(sampl_vec))/fs;
emg_quant = BaseLine2(emg_avg3, 200, 200, 'mn');

end