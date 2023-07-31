function [onset_sampl, onset_time] = emg_onset(emg, fs, W, varargin)
% 
% [onset_sampl, onset_time] = emg_onset(emg, fs, W)
% [onset_sampl, onset_time] = emg_onset(emg, fs, W, th_coeff, Trl)
% 
% *************************************************************************
% * EMG Onset Detection Function                                          *
% *************************************************************************
% 
% Usage:   [onset_sampl, onset_time] = emg_onset(emg, fs, W)
%          [onset_sampl, onset_time] = emg_onset(emg, fs, W, th_coeff, Trl)
% Inputs:
%          'emg': the emg signal
%          'fs': EMG signal's sampling frequency
%          'W': window length for STD calculation
%    (opt) 'th_coeff': threshold coeff for onset detection (default: 1*std)
%    (opt) 'Trl': current trial number
% Outputs:
%          'onset_sampl': sample number of the movement onset
%          'onset_time': corresponding time of movement onset
% Note:
%          an empty bracket [] Must be assigned to not-specified values
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
    Trl = 1;
    th_coeff = 1;
elseif nargin > 3
    if size(varargin, 2) ~= 2
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            Trl = 1;
        else
            if isscalar(varargin{1})
                Trl = varargin{1};
            else
                error('***current trial number has to be a scalar***')
            end
        end
        if isempty(varargin{2})
            th_coeff = 1;
        else
            if isscalar(varargin{2})
                th_coeff = varargin{2};
            else
                error('*** threshold coefficient has to be a scalar***')
            end
        end
    end
end
if isscalar(emg)
    error('***input EMG has to be a double vector or matrix***')
end
if (~(isscalar(fs) && isscalar(W)))
    error('*** sampling frequency and window length have to scalars***')
end
    
%% Detection of an event
for i=1:W:length(emg)-W
    emg_std(((i-1)/W)+1) = std(emg(i:i+W)); % variable size; do not preallocate
end

%% ECG damping
LL = [50, 50];
emg_std_trend = BaseLine2(emg_std, LL(1), LL(2), 'mn');

%% Onset calculation
th = th_coeff*std(emg);
onset_sampl = W*(find(emg_std_trend>=th, 1 ));
onset_time = (onset_sampl)/fs;
disp('*******************************************************************************************************')
fprintf('**  Detected Movement Onset for Trial %d is on sample (%d), corresponding to time (%.2f) seconds  **\n', Trl, onset_sampl, onset_time)
disp('*******************************************************************************************************')

end