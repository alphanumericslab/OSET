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
%        Usage:
%               PLV = PLV_PhaseSeq(phase_sig)
%               PLV = PLV_PhaseSeq(phase_sig1, phase_sig2, phase_sig3, ...)
%        *NOTE: While using the first case, the phase_sig have to be matrix
%               with at least two rows where each row represents a phase
%               signal. In second case, each of the phase_sig1...phase_sign
%               are row vectors of phase sequences.
% 
%             Inputs:
%                       phase_sig(1...n) = phase vectores or a phase matrix
% 
%             Output: 
%                       PLV = PLV matrix
% 
% This program is provided by ESMAEIL SERAJ (esmaeil.seraj09@gmail.com).
% Please make sure to reference BOTH the original studies [1-2] and the 
% OSET [3] to help others find these items.
% 
%     [1] Esmaeil Seraj, Reza Sameni. ”Robust Electroencephalogram Phase 
%         Estimation with Applications in Brain-computer Interface Systems” 
%         Physiological Measurements (2017)
%     [2] Reza Sameni and Esmaeil Seraj, “A Robust Statistical Framework 
%         for Instantaneous Electroencephalogram Phase and Frequency 
%         Analysis” Physiological Measurements (2017)                     
%     [3] R. Sameni, The Open-Source Electrophysiological Toolbox (OSET), 
%         version 3.1 (2014). URL http://www.oset.ir
%         Released under the GNU General Public License
%         Copyright (C) 2012  Reza Sameni
%         Shiraz University, Shiraz, Iran
%         reza.sameni@gmail.com 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.
% 

%%-Checking inputs and assigning default values-%%
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