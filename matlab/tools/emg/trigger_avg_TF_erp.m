function [erp_tf, synch_eeg, trigger_time_sec, time_vec, freq_vec, varargout] = trigger_avg_TF_erp(eeg, fs, onset_time, varargin)
% 
% [erp_tf, synch_eeg, trigger_time_sec, time_vec, freq_vec] = trigger_avg_TF_erp(eeg, fs, onset_time)
% [erp_tf, synch_eeg, trigger_time_sec, time_vec, freq_vec, synch_emg] = trigger_avg_TF_erp(eeg, fs, onset_time, duration, method, emg)
% 
% *************************************************************************
% * Trigger-averaged ERP Time/Frequency Representation                    *
% *************************************************************************
% 
% Usage:    [erp_tf, synch_eeg, trigger_time_sec, time_vec, freq_vec] = trigger_avg_TF_erp(eeg, fs, onset_time)
%           [erp_tf, synch_eeg, trigger_time_sec, time_vec, freq_vec, synch_emg] = trigger_avg_TF_erp(eeg, fs, onset_time, duration, method, emg)
% inputs:
%           'eeg': cell array containing eeg channels of interest from all 
%                  trials
%           'fs': sampling frequency (Hz)
%           'onset_time': vector of onset times where each element is the
%                         onset for that corresponding trial (Seconds)
%     (opt) 'duration': required signal duration after movement onset 
%                       (default: duration = 2 Seconds)
%     (opt) 'method': method for T/F representation, (options: 'STFT', 'CWT'
%                     'NBCH') (default: method = 'STFT')
%     (opt) 'emg': cell array containing emg channel of interest from all 
%                  trials (default: emg = {})
% outputs:
%           'erp_tf': estimated ERP time-frequency map
%           'synch_eeg': synchronized eeg signals based on trigger time
%           'trigger_time_sec': trigger onset flag (Seconds)
%           'time_vec': time vector required for ERP plots
%           'freq_vec': frequency vector required for ERP map plots
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
    method = 'STFT';
    emg = {};
elseif nargin > 3
    if size(varargin, 2) ~= 3
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            duration = 2;
        else
            duration = varargin{1};
        end
        if isempty(varargin{2})
            method = 'STFT';
        else
            method = varargin{2};
        end
        if isempty(varargin{3})
            emg = {};
        else
            emg = varargin{3};
        end
    end
end
if (~(iscell(eeg) && iscell(emg)))
    error('***input EEG/EMG data have to be stored in a cell with each array correspond to one trial***')
end
if (~(isscalar(fs) && isscalar(duration) && isscalar(onset_time)))
    error('***specified sampling frequency and duration have to be scalars***')
end
if ~ischar(method)
    error('***specified time-frequency representation method have to be a string; avaiable options: <STFT> <CWT> <NBCH>***')
else
    if(~(strcmp(freq_band, 'STFT') || strcmp(freq_band, 'CWT') || strcmp(freq_band, 'NBCH')))
        error('***typo in your specified time-frequency representation method string***')
    end
end

%% synchronizing eeg signals based on trigger-time
if isempty(emg)
    [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec] = trigger_synch(eeg, fs, onset_time, duration, emg);
else
    [ensemble_eeg, synch_eeg, trigger_time_sec, time_vec, synch_emg] = trigger_synch(eeg, fs, onset_time, duration, emg);
    varargout = synch_emg;
end
[M, N] = size(ensemble_eeg);

%% calculating the T/F representation using three different methods
if (strcmp(method, 'STFT'))
    % T/F using STFT method
    window = hamming(1024);     % time resolution: about 250ms
    noverlap = 1020;
    f = 2048;                   % frequency resolution
    erp_tf_all = cell(1, M);
    for i=1:M
        [erp_tf_all{i}, freq_vec, time_vec] = spectrogram(ensemble_eeg(i, :), window, noverlap, f, fs);
    end
    erp_tf = zeros(size(erp_tf_all{1}));
    for j=1:M
        erp_tf = erp_tf + erp_tf_all{j};
    end
    erp_tf = erp_tf/M;
    
elseif (strcmp(method, 'CWT'))
    % THIS SECTION IS NOT TESTED
    % T/F using CWT method
    erp_tf_all = cell(1, M);
    scales = 1:64;
    samplingperiod = 1/fs;
    wname = 'db45';  % available options (other versions of MATLAB): 'morse', 'amor', and 'bump'
    for i=1:M
        [erp_tf_all{i}, ~, freq_vec] = cwt(ensemble_eeg(i, :), scales, wname, samplingperiod, 'scal');
    end
    erp_tf = zeros(size(erp_tf_all{1}));
    for j=1:M
        erp_tf = erp_tf + erp_tf_all{j};
    end
    erp_tf = erp_tf/M;
    
elseif (strcmp(method, 'NBCH'))
    % T/F using NBCH method
    f0 = 1:40;          % frequency range of interest
    bw = 4;             % filter bandwidth
    CIC_ord = 10;       % CIC filter order
    bp_filtered_ens = cell(1, M);
    for k=1:M
        for i=1:f0(end)
            bp_filtered_ens{k}(i, :) = BPFilter5(ensemble_eeg(k, :) , f0(i)/fs, bw/fs, CIC_ord);
        end
        bp_filtered_ens{k} = bp_filtered_ens{k}.^2;
    end
    erp_tf_pow = zeros(f0(end), N);
    for j=1:f0(end)
        erp_pow_sum = zeros(1, N);
        for c=1:M
            erp_pow_sum = erp_pow_sum + bp_filtered_ens{c}(j, :);
           erp_tf_pow(j, :) = erp_pow_sum;
        end
    end
    erp_tf_pow = erp_tf_pow/M;
    
    approach = 'mn';
    L = [500 500]; % about 250ms time windows
    erp_tf = BaseLine2(erp_tf_pow, L(1), L(2), approach);
    freq_vec = f0; % frequency vector
    
end

end