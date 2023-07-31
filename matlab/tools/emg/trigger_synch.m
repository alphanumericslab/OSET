function [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec, varargout] = trigger_synch(eeg, fs, onset_time, varargin)
% 
% [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec] = trigger_synch(eeg, fs, onset_time)
% [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec, synch_emg] = trigger_synch(eeg, fs, onset_time, duration, emg)
% 
% *************************************************************************
% * Synchronizing EEG/EMG Signals According to Movement Onset Time        *
% *************************************************************************
% 
% Usage:    [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec] = trigger_synch(eeg, fs, onset_time)
%           [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec, synch_emg] = trigger_synch(eeg, fs, onset_time, duration, emg)
% inputs:
%           'eeg': cell array containing eeg channels of interest from all 
%                  trials (output of 'bdf2mat_main.m' function)
%           'fs': sampling frequency
%           'onset_time': vector of onset times where each element is the
%                         onset for that corresponding trial
%     (opt) 'duration': required signal duration after movement onset
%                       (default: duration = 2 Seconds)
%     (opt) 'emg': cell array containing emg channel of interest from all 
%                  trials (default: emg = {})
% outputs:
%           'ensemble_eeg': ensembles of EEG trials in a matrix (to be used
%                           in 'trigger_avg_erp.m' function)
%           'synch_eeg': synchronized eeg signals based on trigger time
%           'trigger_time_sec': trigger onset flag (seconds)
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
if nargin < 3
    error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 3
    duration = 2;
    emg = {};
elseif nargin > 3
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

%% eeg synchronization
onset_sampl = onset_time*fs;
synch_eeg = cell(1, length(eeg));
trigger_time = min(onset_sampl);
trigger_time_sec = trigger_time/fs;
post_move_time = duration*fs;
for i=1:length(eeg)
    sampl_vec = onset_sampl(i)-trigger_time+1:onset_sampl(i)+post_move_time;
    synch_eeg{i} = eeg{i}(sampl_vec);
end
time_vec = (sampl_vec)/fs;
ensemble_eeg = zeros(length(synch_eeg), length(sampl_vec));
for j=1:length(synch_eeg)
    ensemble_eeg(j, :) = synch_eeg{i};
end

%% emg synchronization
if ~isempty(emg)
    synch_emg = cell(1, length(emg));
    for i=1:length(emg)
        sampl_vec = onset_sampl(i)-trigger_time+1:onset_sampl(i)+post_move_time;
        synch_emg{i} = emg{i}(sampl_vec);
    end
    varargout = synch_emg;
end

end