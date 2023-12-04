function pwplv = PWPLV(eeg, fs, onset_time, varargin)
% 
% pwplv = PWPLV(eeg, fs, onset_time)
% pwplv = PWPLV(eeg, fs, onset_time, freq_rng, duration, pertnum, plot_flag)
% 
% *************************************************************************
% * Pair-wise PLV dynamics estimated within 1sec time-steps for any       *
% * arbitrary time range and all electrode pairs                          *
% *************************************************************************
% 
% Usage:    pwplv = PWPLV(eeg, fs, onset_time)
%           pwplv = PWPLV(eeg, fs, onset_time, freq_rng, duration, pertnum, plot_flag)
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
%     (opt) 'pertnum': number of perturbations while using the TFS phase
%                      estimation method (default: pertnum = 100)
%     (opt) 'plot_flag': decide to visualize the results or not. options: 
%                        'plot', 'noplot' (default: plot_flag = 'plot')
% outputs:
%           'pwplv': estimated pairwise phase locking value (PLV) between
%                    all possible electrode pairs
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
    pertnum = 100;
    plot_flag = 'plot';
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
            pertnum = 100;
        else
            pertnum = varargin{3};
        end
        if isempty(varargin{4})
            plot_flag = 'plot';
        else
            plot_flag = varargin{4};
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
if ischar(plot_flag)
    if (~(strcmp(plot_flag, 'plot') || strcmp(plot_flag, 'noplot')))
		error('***typo in your specified plot_flag string***')
    end
else
    error('***plot_flag has to be a string***')
end

%% initialization and parameter specification
f0 = (freq_rng(2)+freq_rng(1))/2;             % center frequency for filter
bw = freq_rng(2)-freq_rng(1);                         % frequency bandwidth
time_steps = (duration(1):duration(2)).*fs;   % time steps for MSC measures

%% pairwise phase locking value (PLV) estimation within time steps
% check if parallel pool has been created
poolsize = 4; % if you have more cores available GOOD FOR YOU, change it here
p = gcp('nocreate');                    % If no pool, do not create new one
if isempty(p)
    parpool(poolsize)
end

onset_sampl = onset_time*fs;
plv_gen = cell(1, length(eeg));
for i=1:length(eeg)
    bp_eeg = eeg{i};
    [m, ~] = size(bp_eeg);
    phase_eeg = zeros(size(bp_eeg));
    
    % phase estimation by TFP method
    parfor ii=1:m
        [phase_eeg(ii, :), ~, ~, ~] = phase_est(bp_eeg(ii, :), fs, f0, bw, pertnum);
    end
    
    % PLV estimation
    PLV_all = cell(1, length(time_steps)-1);
    for k=1:length(time_steps)-1
        time_vec = onset_sampl+time_steps(k):onset_sampl+time_steps(k+1);
        plv = zeros(m);
        for j=1:m
            for c=1:m
                plv(j, c) = PLV_PhaseSeq(phase_eeg(j, time_vec), phase_eeg(c, time_vec));
            end
        end
        PLV_all{k} = plv;
    end
    plv_gen{i} = PLV_all;
end

% averaging PLV measures across trials
pwplv = cell(1, length(time_steps)-1);
for i=1:length(time_steps)-1
    sum_plv = zeros(m);
    for j=1:length(eeg)
        sum_plv = sum_plv + plv_gen{j}{i};
    end
    mean_plv = sum_plv/length(eeg);
    pwplv{i} = mean_plv;
end

%% visualizing the results
if (strcmp('plot', plot_flag))
    figure
    L = length(pwplv);
    switch L
        case {1, 2, 3}
            ind1 = 1;
            ind2 = L;
        case {4}
            ind1 = 2;
            ind2 = 2;
        case {5, 6}
            ind1 = 2;
            ind2 = 3;
        case {7, 8, 9}
            ind1 = 3;
            ind2 = 3;
    end
    for i=1:length(pwplv)
        subplot(ind1, ind2, i)
        contourf(pwplv{i}, 8)
        colormap jet
        xlabel(sprintf('(%d) to (%d) Sec', time_steps(i)/fs, time_steps(i+1)/fs));
        ylabel(sprintf('Pairwise PLV (%d-%dHz)', freq_rng(1), freq_rng(2)))
    end
elseif (strcmp('noplot', plot_flag))
    warning('Result visualization is OFF, change the setting if you wish to visualize the connectivity results!!')
end

end