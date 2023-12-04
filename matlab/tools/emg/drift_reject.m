function sig = drift_reject(raw_sig, L1, varargin)
% 
% sig = drift_reject(raw_sig, L1)
% sig = drift_reject(raw_sig, L1, L2, approach)
% 
% *************************************************************************
% * Drift (Baseline Wander) Cancellation from Biological Recordings       *
% *************************************************************************
% 
% Usage:    sig = drift_reject(raw_sig, L1)
%           sig = drift_reject(raw_sig, L1, L2, approach)
% inputs:
%           'raw_sig': matrix or vector of raw recordings
%           'L1': first stage window length
%     (opt) 'L2': second stage window length (default: L2 = L1)
%     (opt) 'approach': filtering approach, available options: 'md', 'mn'
%                       (default: approach = 'mn')
% outputs:
%           'sig': matrix or vector of drift-rejected signals
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
if nargin < 2
	error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 2
    L2 = L1;
    approach = 'mn';
elseif nargin > 2
    if (size(varargin, 2) ~= 2)
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin{1})
            L2 = L1;
        else
            if isscalar(varargin{1})
                L2 = varargin{1};
            else
                error('*** L2 has to be a scalar***')
            end
        end
        if isempty(varargin{2})
            approach = 'mn';
        else
            if ischar(varargin{2})
                approach = varargin{2};
            else
                error('***approach name has to be a character***')
            end
        end
    end
end
if (isscalar(raw_sig))
	error('***input raw signals needs to be a double vector or matrix***')
end
if (~(isscalar(L1) && isscalar(L2)))
	error('***window lengthes have to be scalars***')
end
if ischar(approach)
	if (~(strcmp(approach, 'mn') || strcmp(approach, 'md')))
		error('***typo in your specified filtering approach string***')
	end
else
	error('***specified approach name has to be a string***')
end

%% drift rejection
drift_baseline = BaseLine2(raw_sig, L1, L2, approach);
sig = raw_sig - drift_baseline;

end