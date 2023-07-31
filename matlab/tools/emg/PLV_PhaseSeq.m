function PLV = PLV_PhaseSeq(varargin)
% 
% PLV = PLV_PhaseSeq(phase_sig)
% PLV = PLV_PhaseSeq(phase_sig1, phase_sig2, phase_sig3, ...)
% 
% *************************************************************************
% * Calculating Phase Locking Value (PLV) matrix (Pairwise) using phase   *
% * sequences. Refer to the User Guide for further details.               *
% *************************************************************************
% 
% Usage:    PLV = PLV_PhaseSeq(phase_sig)
%           PLV = PLV_PhaseSeq(phase_sig1, phase_sig2, phase_sig3, ...)
% inputs:
%           phase_sig(1...n) = phase vectores or a phase matrix
% outputs:
%           PLV = PLV matrix
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
if nargin<1
    error('***Calculating Pairwise PLV Requires at Least TWO Phase Sequences***')
end
if nargin<2
    [aa, bb] = size(varargin{1});
    if aa==1
        error('***Calculating Pairwise PLV Requires at Least TWO Phase Sequences***')
    else
        L = bb;
        phase_sig = cell(1, aa);
        for j=1:aa
            phase_sig{j} = varargin{1}(j, :);
        end
    end
else
    L = length(varargin{1});
    phase_sig = cell(1, length(varargin));
    for i=1:length(varargin)
        phase_sig{i} = varargin{i};
    end
end

%% PLV estimation
[~, n1] = size(phase_sig);
phase_sig = repmat(phase_sig, n1, 1);
[~, nn] = size(phase_sig);

phase_diff = cell(nn);
for k1=1:nn
    for k2=1:nn
        phase_diff{k1, k2} = phase_sig{k1, k1} - phase_sig{k2, k2};
    end
end

PLV = zeros(nn);
for l1=1:nn
    for l2=1:nn
        PLV(l1, l2) = abs(sum(exp(1j*(phase_diff{l1, l2}))))/L;
    end
end
[PLVm, ~] = size(PLV);
if PLVm==2
	PLV = PLV(1, 2);
end

end