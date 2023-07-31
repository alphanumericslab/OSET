function tcplv = TCPLV(eeg, fs, onset_time, varargin)
% 
% tcplv = TCPLV(eeg, fs, onset_time)
% tcplv = TCPLV(eeg, fs, onset_time, freq_rng, duration, pairofint, pertnum)
% 
% *************************************************************************
% * PLV temporal dynamics estimated within 1sec time-steps for any        *
% * arbitrary time range before and after movement onset                  *
% *************************************************************************
% 
% Usage:    tcplv = TCPLV(eeg, fs, onset_time)
%           tcplv = TCPLV(eeg, fs, onset_time, freq_rng, duration, pairofint, pertnum)
% inputs:
%           'eeg': cell array containing eeg channels of interest from all 
%                  trials
%           'fs': sampling frequency (Hz)
%           'onset_time': vector of onset times (Seconds)
%     (opt) 'freq_rng': [a, b] form double vector where 'a' and 'b' are 
%                       edges of the frequency band of interest (default:
%                       freq_rng = [12, 32])
%     (opt) 'duration': [-a, b] form double vector where 'a' is time required 
%                       duration prior to movement onset and 'b' is the required 
%                       time duration after the movement onset in seconds (
%                       default: duration = [-3, 2])
%     (opt) 'pairofint': channel pair of interest for PLV time-course (
%                        default: pairofint = 'all')
%     (opt) 'pertnum': number of perturbations while using the TFS phase
%                      estimation method (default: pertnum = 100)
% outputs:
%           'tcplv': estimated time-course phase locking value (PLV) dynamics
%                    between eeg pairs of interest
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
    freq_rng = [12, 32];
    duration = [-3, 2];
    pairofint = 'all';
    pertnum = 100;
elseif nargin > 3
    if size(varargin, 2) ~= 4
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            freq_rng = [12, 32];
        else
            freq_rng = varargin{1};
        end
        if isempty(varargin{2})
            duration = [-3, 2];
        else
            duration = varargin{2};
        end
        if isempty(varargin{3})
            pairofint = 'all';
        else
            pairofint = varargin{3};
        end
        if isempty(varargin{4})
            pertnum = 100;
        else
            pertnum = varargin{4};
        end
    end
end
if (~(isscalar(fs) && isscalar(onset_time) && isscalar(pertnum)))
    error('***fs, trigger time and number of perturbations have to be scalars***')
end
if(isscalar(freq_rng) || isscalar(duration))
    error('***frequency range and duration of interest have to be double vectors in form of [a, b]. Refer to manual for details***')
end
if ~iscell(eeg)
    error('***input EEG signal has to be stored in a cell array. Refer to manual for details***')
end
if (~(isa(pairofint, 'double') || ischar(pairofint)))
    error('***wrong format for pariofint variable***')
else
    if isscalar(pairofint)
        error('***wrong format for pariofint variable***')
    end
end

%% initialization and parameter specification
f0 = (freq_rng(2)+freq_rng(1))/2;             % center frequency for filter
bw = freq_rng(2)-freq_rng(1);                 % frequency bandwidth
time_steps = (duration(1):duration(2)).*fs;   % time steps for MSC measures (250ms)

%% pairwise phase locking value (PLV) estimation within time steps
% check if parallel pool has been created
poolsize = 4; % if you have more cores available GOOD FOR YOU, change it here
p = gcp('nocreate');                    % If no pool, do not create new one
if isempty(p)
    parpool(poolsize)
end

onset_sampl = onset_time*fs;
if ischar(pairofint) % in case you want a colored map of tcplv between C3 and every other electrodes
    plv_gen = cell(length(eeg));
    
    for i=1:length(eeg)
        % phase estimation by TFP method
        sig_eeg = eeg{i};
        phase_eeg = zeros(size(sig_eeg));
        parfor ii=1:size(sig_eeg, 1)
            [phase_eeg(ii, :), ~, ~, ~] = phase_est(sig_eeg(ii, :), fs, f0, bw, pertnum);
        end
        
        plv = zeros(size(phase_eeg, 1), length(time_steps)-1);
        for j=1:size(phase_eeg, 1)
            % PLV estimation
            for k=1:length(time_steps)-1
                time_vec = onset_sampl+time_steps(k):onset_sampl+time_steps(k+1);
                plv(j, k) = PLV_PhaseSeq(phase_eeg(6, time_vec), phase_eeg(j, time_vec)); % 6: C3
            end
        end
        plv_gen{i} = plv;
    end

    % averaging PLV measures across trials
    plv_sum = zeros(size(plv_gen{1}));
    for ll=1:length(plv_gen)
        plv_sum = plv_sum + plv_gen{ll};
    end
    tcplv = plv_sum/length(plv_gen);
    
elseif isdouble(pairofint) % in case you only want a single curve between two electrodes
    
    plv_gen = zeros(length(eeg), length(time_steps)-1);
    for i=1:length(eeg)
        bp_eeg = [eeg{i}(pairofint(1), :); eeg{i}(pairofint(2), :)];
        phase_eeg = zeros(size(bp_eeg));
        % phase estimation by TFP method
        parfor ii=1:size(bp_eeg, 1)
            [phase_eeg(ii, :), ~, ~, ~] = phase_est(bp_eeg(ii, :), fs, f0, bw, pertnum);
        end

        % PLV estimation
        plv = zeros(1, length(time_steps)-1);
        for k=1:length(time_steps)-1
            time_vec = onset_sampl+time_steps(k):onset_sampl+time_steps(k+1);
            plv(k) = PLV_PhaseSeq(phase_eeg(1, time_vec), phase_eeg(2, time_vec));
        end
        plv_gen(i, :) = plv;
    end

    % averaging PLV measures across trials
    tcplv = mean(plv_gen);
    
else
    error('***The entered pair-of-interest format is not correct***')
end

end