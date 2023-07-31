function [PR, PS, PL, PDV, PD] = Phase_Features_MultiCh(meth, sig1, sig2, fs, f0, varargin)
% 
% [PR, PS, PL, PDV, PD] = Phase_Features_MultiCh(meth, sig1, sig2, fs, f0)
% [PR, PS, PL, PDV, PD] = Phase_Features_MultiCh(meth, sig1, sig2, fs, f0, WS, ndft, pertnum)
% 
% *************************************************************************
% * Calculating popular phase related quantities, i.e. Phase Shift (PS),  *
% * Phase Lock (PL), Phase Reset (PR), Phase Difference (PD) and Phase    *
% * Difference Derivatives (PDV) in MULTI-Channel mode. Refer to the User *
% * Guide for further details.                                            *
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
%                       sig1 : input signal #1
%                       sig2 : input signal #2
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
%                       PS : Phase Shift Events
%                       PL : Phase Lock Events
%                       PR : Phase Resetting Events
%                       PD : Phase Difference
%                       PDV : Phase Difference Derivatives
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
sig1 = sig1(:)';
sig2 = sig2(:)';
if nargin<5
    error('***Not Enough Input Arguments: Please Check the Instructions***')
elseif nargin==5
    WS_trad = 1;
    WS_pert = 1;
    ndft = 100;
    pertnum = 100;
elseif nargin>5
    if nargin>8
        error('***Too Many Input Arguments: Please Check the Instructions***')
    end
    if nargin~=8
        error('***An Empty Bracket [] Must be Used for the Values NOT Being Specified***')
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
end

[m1, n1] = size(sig1);
[m2, n2] = size(sig2);
if m1~=m2 || n1~=n2
    error('***Input Signals Dimension Mismatch***')
end
if m1~=1 || m2~=1
    error('***Input Signals Must Contain Only One Channel of Information***')
end
meth1 = 'Trad';
meth2 = 'ZPPP';

%%-Phase Calculation-%%
if (strcmp(meth1, meth)==1) 
    [phase_sig1, ~, ~] = Phase_Ex_Trad(sig1, fs, WS_trad, f0, ndft);
    [phase_sig2, ~, ~] = Phase_Ex_Trad(sig2, fs, WS_trad, f0, ndft);
    flg = 0;

elseif (strcmp(meth2, meth)==1)
    [phase_sig1, ~, ~] = Phase_Ex_ZPPP(sig1, fs, WS_pert, f0, ndft, pertnum);
    [phase_sig2, ~, ~] = Phase_Ex_ZPPP(sig2, fs, WS_pert, f0, ndft, pertnum);
    flg = 1;
    
else
    error('***Please Enter the Name of the Phase Estimation Method Correctly (Available methods are "ZPPP" and "Trad")***')
end
phase_sig1 = phase_sig1(f0, :);
phase_sig2 = phase_sig2(f0, :);

%%-Phase Features Estimation-%%
PD = phase_sig1 - phase_sig2;
PDV = [diff(PD), 0];

if flg == 0
    Th_PS = 0.2;    % Threshold for PS and PL events (Radians): Subjective
else
    Th_PS = 0.02;   % Threshold for PS and PL events (Radians): Subjective
end

PS = zeros(m1 ,n1);
for j=1:n1-1
    if PDV(j)>=Th_PS || PDV(j)<=-Th_PS
        PS(j) = 1+Th_PS/100;
    end
end
PL = ~PS;
PR = ones(m1, n1)+Th_PS/100;        % shifting a bit to have a better view
ind = find(PL~=0);
for c=1:length(ind)-1
    if ind(c+1)~=ind(c)+1
        PR(1, ind(c)+1:end) = -PR(1, ind(c)+1:end);
    end
end
PL = PL+Th_PS/100;

end
