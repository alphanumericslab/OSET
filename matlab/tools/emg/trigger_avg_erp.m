function [erp, synch_eeg, trigger_time_sec, time_vec, varargout] = trigger_avg_erp(eeg, fs, freq_band, onset_time, varargin)
% 
% [erp, synch_eeg, trigger_time_sec, time_vec] = trigger_avg_erp(eeg, fs, freq_band, onset_time)
% [erp, synch_eeg, trigger_time_sec, time_vec, synch_emg] = trigger_avg_erp(eeg, fs, freq_band, onset_time, duration, emg)
% 
% *************************************************************************
% * Trigger-averaged ERP time-course estimation function                  *
% *************************************************************************
% 
% Usage:    [erp, synch_eeg, trigger_time_sec, time_vec] = trigger_avg_erp(eeg, fs, freq_band, onset_time)
%           [erp, synch_eeg, trigger_time_sec, time_vec, synch_emg] = trigger_avg_erp(eeg, fs, freq_band, onset_time, duration, emg)
% inputs:
%           'eeg': cell array containing eeg channels of interest from all 
%                  trials
%           'fs': sampling frequency (Hz)
%           'freq_band': frequency band of interest, options: 'delta', 
%                        'theta', 'alpha', 'beta' or 'gamma'
%           'onset_time': vector of onset times where each element is the
%                         onset for that corresponding trial (Seconds)
%     (opt) 'duration': required signal duration after movement onset 
%                       (default: duration = 2 Seconds)
%     (opt) 'emg': cell array containing emg channel of interest from all 
%                  trials (default: emg = {})
% outputs:
%           'erp': extracted ERP signal
%           'synch_eeg': synchronized eeg signals based on trigger time
%           'trigger_time_sec': trigger onset flag (Seconds)
%           'time_vec': time vector required for ERP plots
%           'synch_emg': synchronized emg signals based on trigger time
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
if nargin < 4
    error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 4
    duration = 2;
    emg = {};
elseif nargin > 4
    if size(varargin, 2) ~= 2
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            duration = 2;
        else
            duration = varargin{1};
        end
        if isempty(varargin{2})
            emg = {};
        else
            emg = varargin{2};
        end
    end
end
if (~(iscell(eeg) && iscell(emg)))
    error('***input EEG/EMG data have to be stored in a cell with each array correspond to one trial***')
end
if (~(isscalar(fs) && isscalar(duration) && isscalar(onset_time)))
    error('***specified sampling frequency and duration have to be scalars***')
end
if ~ischar(freq_band)
    error('***specified frequency band have to be a string; avaiable options: <delta> <theta> <alpha> <beta> <gamma>***')
else
    if(~(strcmp(freq_band, 'delta') || strcmp(freq_band, 'theta') || strcmp(freq_band, 'alpha') || strcmp(freq_band, 'beta') || strcmp(freq_band, 'gamma')))
        error('***typo in your specified frequency band string***')
    end
end

%% synchronizing eeg signals based on trigger-time
if isempty(emg)
    [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec] = trigger_synch(eeg, fs, onset_time, duration, emg);
else
    [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec, synch_emg] = trigger_synch(eeg, fs, onset_time, duration, emg);
    varargout = synch_emg;
end

%% trigger-averaged ERP estimation
if (strcmp(freq_band, 'delta')) % frequency range of interest for "Delta" band: 1-4Hz
    f0 = 2.5; bw = 3; 
elseif (strcmp(freq_band, 'theta')) % frequency range of interest for "Theta" band: 4-8Hz
    f0 = 6; bw = 4;
elseif (strcmp(freq_band, 'alpha')) % frequency range of interest for "Alpha" band: 8-12Hz
    f0 = 10; bw = 4;
elseif (strcmp(freq_band, 'beta')) % frequency range of interest for "Beta" band: 12-32Hz
    f0 = 22; bw = 20;
elseif (strcmp(freq_band, 'gamma')) % frequency range of interest for "Gamma" band: 32-80Hz
    f0 = 56; bw = 48;
else
	error('The Name of the EEG Frequency Band Misspelled, Check for Typos..!')
end
CIC_ord = 5; % you can change the filter order here
bp_filtered_eeg = BPFilter5(ensemble_eeg , f0/fs, bw/fs, CIC_ord);
erp_env = mean(bp_filtered_eeg.^2);

% [erp, erp_vec_plot] = sig_trend(erp_env);
% erp_vec_plot = erp_vec_plot/fs;
approach = 'mn';
L = [500 500];          % about 250ms time windows
erp = BaseLine2(erp_env, L(1), L(2), approach);

end