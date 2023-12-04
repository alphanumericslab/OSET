function PLV = PLV_RawSig(meth, sig1, fs, f0, varargin)
% 
% PLV = PLV_RawSig(meth, sig1, fs, f0)
% PLV = PLV_RawSig(meth, sig1, fs, f0, WS, ndft, pertnum, sig2)
% 
% *************************************************************************
% * Calculating Phase Locking Value (PLV) matrix (Pairwise) using raw     *
% * signals. Refer to the User Guide for further details.                 *
% *************************************************************************
% 
%             Usage:
%                       meth = 'Trad' *** calculating the phase features
%                                         through using the conventional 
%                                         approach for phase estimation
%                                         based on analytic representation.
%                       meth = 'ZPPP' *** calculating the phase features
%                                         through using the Zero-pole 
%                                         perturbation approach for phase 
%                                         estimation [1].
% 
%             Inputs:
%                       meth : the utilized method for phase estimation
%                       sig1 : input signal #1 *** Could be a matrix (while
%                              calculating PLV in multi-channel signals) or
%                              a vector (while calculating PLV between two
%                              separate signals).
%            (optional) sig2 : input signal #2 *** This option is provided
%                              in case that someone needs to calculate the
%                              PLV between two separate signals.
%                       fs : sampling frequency (Hz)
%                       f0 : the frequency of interest (Hz)
%            (optional) WS : Stop-band frequency (Hz) >> (BW = WS/2)
%            (optional) ndft: number of frequency bins
%            (optional) pertnum: number of attempts for perturbing filter's
%                                zeros and poles *** Not required while
%                                using 'Trad' method for phase estimation
% 
%             Defaults:
%                       WS = 1(Hz)
%                       ndft = 100
%                       pertnum = 100
% 
%             *NOTE: While specifying a value to one of the parameters
%                    having default values, an empty bracket [] must be
%                    used for non-specified parameters. If you're not using
%                    these parameters, empty bracket is not required.
%             Outputs: 
%                       PLV : PLV matrix 
% 
%             *NOTE: While trying to calculate the PLV between two separate
%                    signals, if you don't need a PLV matrix, just simply
%                    use one of the non-diagonal values in estimated PLV
%                    matrix as the PLV between signals.
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
if nargin<4
    error('***Not Enough Input Arguments: Please Check the Instructions***')
elseif nargin==4
    WS_trad = 1;
    WS_pert = 1;
    ndft = 100;
    pertnum = 100;
elseif nargin>4
        if nargin>8
            error('***Too Many Input Arguments: Please Check the Instructions***')
        end
        if nargin<7
            error('***An Empty Bracket [] Must be Used for the Values NOT Being Specified***')
        end
end
if nargin==7
    [m1, L] = size(sig1);
    if isempty(varargin{1})
        WS_trad = 1;
        WS_pert = 1;
    else
        WS_trad = varargin{1};
        WS_pert = varargin{1};
    end
    if isempty(varargin{2})
        ndft = 100;
    else
        ndft = varargin{2};
    end
    if isempty(varargin{3})
        pertnum = 100;
    else
        pertnum = varargin{3};
    end
elseif nargin==8
    sig2 = varargin{4};
    sig2 = sig2(:)';
    [m2, n2] = size(sig2);
    [m1, n1] = size(sig1);
    if m1~=1 || m2~=1
        error('***In Multi-signal Case, Signals Must Contain Only One Channel of Information***')
    end
    if m1~=m2 || n1~=n2
        error('***Input Signals Dimension Mismatch***')
    end
	if isempty(varargin{1})
        WS_trad = 1;
        WS_pert = 1;
    else
        WS_trad = varargin{1};
        WS_pert = varargin{1};
    end
    if isempty(varargin{2})
        ndft = 100;
    else
        ndft = varargin{2};
    end
    if isempty(varargin{3})
        pertnum = 100;
    else
        pertnum = varargin{3};
    end
    sig1 = [sig1; sig2];
    [m1, L] = size(sig1);
elseif nargin<=7
    [m1, L] = size(sig1);
    if m1==1
        error('***In Multi-channel Case, Input Signal Must Contain at Least TWO Channels of Information***')
    end
end
sig = cell(1, m1);
for i=1:m1
    sig{i} = sig1(i, :);
end
meth1 = 'Trad';
meth2 = 'ZPPP';

%%-phase calculation-%%
phase_sig = cell(1, m1);
if (strcmp(meth1, meth)==1)
    for j=1:m1
        [phase_sig1, ~, ~] = Phase_Ex_Trad(sig{j}, fs, WS_trad, f0, ndft);
        phase_sig{j} = phase_sig1(f0, :);
    end
            
elseif (strcmp(meth2, meth)==1)
    for c=1:m1
        [phase_sig1, ~, ~] = Phase_Ex_ZPPP(sig{c}, fs, WS_pert, f0, ndft, pertnum);
        phase_sig{c} = phase_sig1(f0, :);
    end
    
else
error('***Please Enter the Name of the Phase Estimation Method Correctly (Available methods are "ZPPP" and "Trad")***')
end
phase_sig = repmat(phase_sig, m1, 1);
[~, nn] = size(phase_sig);

%%-calculating the PLV matrix-%%
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
if m1==2
	PLV = PLV(1, 2);
end

end