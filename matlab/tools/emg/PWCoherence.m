function pwcoher = PWCoherence(eeg, fs, onset_time, varargin)
% 
% pwcoher = PWCoherence(eeg, fs, onset_time)
% pwcoher = PWCoherence(eeg, fs, onset_time, freq_rng, duration, plot_flag)
% 
% *************************************************************************
% * Pair-wise MSC dynamics estimated within 1sec time-steps for any       *
% * arbitrary time range and all electrode pairs                          *
% *************************************************************************
% 
% Usage:    pwcoher = PWCoherence(eeg, fs, onset_time)
%           pwcoher = PWCoherence(eeg, fs, onset_time, freq_rng, duration, plot_flag)
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
%     (opt) 'plot_flag': decide to visualize the results or not. options: 
%                        'plot', 'noplot' (default: plot_flag = 'plot')
% outputs:
%           'pwcoher': estimated pairwise magnitude squared coherence (MSC)
%                      between all possible electrode pairs
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
    plot_flag = 'plot';
elseif nargin > 3
    if size(varargin, 2) ~= 3
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
            plot_flag = 'plot';
        else
            plot_flag = varargin{3};
        end
    end
end
if (~(isscalar(fs) && isscalar(onset_time)))
    error('***fs and trigger time have to be scalars. Refer to manual for more details***')
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
window = hamming(1024);                     % window for MSC calculation
nfft = 2048;                                % frequency bins
noverlap = 1;                               % number of overlapping samples
% ch = ceil(nfft.*freq_rng/(fs/2))/2;
f0 = (freq_rng(2)+freq_rng(1))/2;           % center frequency for filter
bw = freq_rng(2)-freq_rng(1);               % bandwidth of frequency filter
CIC_ord = 5;                                % order of frequency filter
time_steps = (duration(1):duration(2)).*fs; % time steps for MSC measures

%% pairwise magnitude squared coherence estimation within time steps
onset_sampl = onset_time*fs;
coher = cell(1, length(eeg));
for i=1:length(eeg)
    % frequency filtering within specified frequency range
    bp_filtered_eeg = BPFilter5(eeg{i} , f0/fs, bw/fs, CIC_ord);
    [m, ~] = size(bp_filtered_eeg);
    
    % coherence estimation
    MSC = cell(1, length(time_steps)-1);
    for k=1:length(time_steps)-1
        time_vec = onset_sampl+time_steps(k):onset_sampl+time_steps(k+1);
        Cxy = zeros(m);
        for j=1:m
            for c=1:m
                Cxy(j, c) = mean(mscohere(bp_filtered_eeg(j, time_vec), bp_filtered_eeg(c, time_vec),...
                    window, noverlap, nfft, fs, 'twosided'));
            end
        end
        MSC{k} = Cxy;
    end
    coher{i} = MSC;
end

%% averaging MSC measures across trials
pwcoher = cell(1, length(time_steps)-1);
for i=1:length(time_steps)-1
    sum_coher = zeros(m);
    for j=1:length(eeg)
        sum_coher = sum_coher + coher{j}{i};
    end
    mean_coher = sum_coher/length(eeg);
    pwcoher{i} = mean_coher;
end

%% visualizing the results
if (strcmp('plot', plot_flag))
    figure
    L = length(pwcoher);
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
    for i=1:length(pwcoher)
        subplot(ind1, ind2, i)
        contourf(pwcoher{i}, 8)
        colormap jet
        xlabel(sprintf('(%d) to (%d) Sec', time_steps(i)/fs, time_steps(i+1)/fs));
        ylabel(sprintf('Pairwise MSC (%d-%dHz)', freq_rng(1), freq_rng(2)))
    end
elseif (strcmp('noplot', plot_flag))
    warning('Result visualization is OFF, change the setting if you wish to visualize the connectivity results!!')
end

end