function [eeg_data, emg_data, fs, varargout] = bdf2mat_main(trl_num, elec_num, emg_ch, eeg_ch, filename, varargin)
% 
% [eeg_data, emg_data, fs] = bdf2mat_main(trl_num, elec_num, emg_ch, eeg_ch, filename)
% [eeg_data, emg_data, fs, emg_onset_sampl, emg_onset_time] = bdf2mat_main(trl_num, elec_num, emg_ch, eeg_ch, filename, drift_flag, onset_flag)
% 
% *************************************************************************
% * Reading EEG/EMG signal from bdf files, preprocessing and conditioning *
% *************************************************************************
% 
% Usage:    [eeg_data, emg_data, fs] = bdf2mat_main(trl_num, elec_num, emg_ch, eeg_ch, filename)
%           [eeg_data, emg_data, fs, emg_onset_sampl, emg_onset_time] = bdf2mat_main(trl_num, elec_num, emg_ch, eeg_ch, filename, drift_flag, onset_flag)
% inputs:
%           'trl_num': number of trials
%           'elec_num': number of electrodes used to record EEG
%           'emg_ch': EMG channel of interest
%           'eeg_ch': scalar or double vector of EEG channels of interest
%           'filename': file name format as a string
%     (opt) 'drift_flag': baseline drift rejection flag, available options:
%                         'drift', 'nodrift'(default: drift_flag = 'drift')
%     (opt) 'onset_flag': onset detection flag, available options: 1, 0
%                         (default = 1)
% outputs:
%           'eeg_data': preprocessed EEGs of the selected channel from all
%                       trials
%           'emg_data': preprocessed EMGs of the selected channel from all
%                       trials
%           'fs': sampling frequency
%           'emg_onset_sampl': sample number of the movement onset
%           'emg_onset_time': corresponding time of movement onset
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
if nargin < 5
    error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 5
    drift_flag = 'drift';
    onset_flag = 1;
elseif nargin > 5
    if size(varargin, 2) ~= 2
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            drift_flag = 'drift';
        else
            drift_flag = varargin{1};
        end
        if isempty(varargin{2})
            onset_flag = 1;
        else
            onset_flag = varargin{2};
        end
    end
end
if (~(isscalar(trl_num) && isscalar(elec_num) && isscalar(onset_flag)))
    error('***trial and electrode numbers have to be a scalar***')
end
if (~(isa(emg_ch, 'double') && isa(eeg_ch, 'double')))
    error('***EMG and EEG channels of interest have to be a doubles (i.e. integer or vector of indices)***')
end
if (~(ischar(filename) && ischar(drift_flag)))
    error('***filename and drift_flag have to be strings***')
else
    if (~(strcmp(drift_flag, 'drift') || strcmp(drift_flag, 'nodrift')))
		error('***typo in your specified drift_flag string***')
    end
end

%% *.bdf to *.mat transformation
emg_data = cell(1, trl_num);
eeg_data = cell(1, trl_num);
emg_onset_sampl = zeros(1, trl_num);
emg_onset_time = zeros(1, trl_num);
L = [700, 500];         % time window needed here has to be about 1-2 sec
approach = 'mn';
for i=1:trl_num
    % reading data and separating EEG & EMG signals
% [data,numChan,labels,txt,fs,gain,prefiltering,ChanDim] = eeg_read_bdf(sprintf(filenameformat, i), 'all', 'n');
    [data, ~, ~, ~, fs, ~, ~, ~] = eeg_read_bdf(sprintf(filename, i), 'all', 'n');
    eeg_all = data(1:elec_num ,:);
    if (isscalar(eeg_ch) || isscalar(emg_ch))
        emg = data(emg_ch, :);
        eeg = eeg_all(eeg_ch, :);
    elseif (isa(eeg_ch, 'double') || isa(emg_ch, 'double'))
        [~, nn] = size(data);
        emg = zeros(length(emg_ch), nn);
        eeg = zeros(length(eeg_ch), nn);
        for lll=1:length(emg_ch)
            emg(lll, :) = data(emg_ch(lll), :);
        end
        for ll=1:length(eeg_ch)
            eeg(ll, :) = eeg_all(eeg_ch(ll), :);
        end
    end
%     ref = mean(data(1:16, :));
%     ref = (data(271, :) + data(272, :))/2;
%     eeg = eeg-ref;
    
    if (strcmp(drift_flag, 'drift'))
        % preprocessing the data: BaseLineWander Removal
        emg = drift_reject(emg, L(1), L(2), approach);
        eeg = drift_reject(eeg, L(1), L(2), approach);
    elseif (strcmp(drift_flag, 'nodrift'))
        emg = drift_reject(emg, L(1), L(2), approach);
        if onset_flag == 1
            window = 30;
            [emg_onset_sampl(i), emg_onset_time(i)] = emg_onset(emg, fs, window);
        end
    else
        error('You need to choose an option to deal with DRIFT')
    end
    
    % EMG onset estimation
    if onset_flag == 1
        window = 30;
        [emg_onset_sampl(i), emg_onset_time(i)] = emg_onset(emg, fs, window, i);
    end
    
    % saving EEG & EMG signal
    emg_data{i} = emg;
    eeg_data{i} = eeg;
    eeg = [];
end
if onset_flag == 1
    varargout{1} = emg_onset_sampl;
    varargout{2} = emg_onset_time;
end

end