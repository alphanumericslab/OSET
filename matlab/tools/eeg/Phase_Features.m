function varargout = Phase_Features(phase_sig1,  fs, Th, varargin)
% 
% varargout = Phase_Features(phase_sig1,  fs, Th)
% varargout = Phase_Features(phase_sig1,  fs, Th, phase_sig2)
% 
% *************************************************************************
% * Calculating popular phase related quantities, i.e. Phase Shift (PS),  *
% * Phase Lock (PL), Phase Reset (PR), Phase Difference (PD) and          *
% * Instantaneous Frequency (IF) in SINGLE or MULTI-Channel modes using   *
% * phase sequences. Refer to the User Guide for further details.         *
% *************************************************************************
% 
% Usage:
%   Multi-channel :  
%   [PR, PS, PL, PDV, PD] = Phase_Features(phase_sig1,  fs, Th, phase_sig2)
%   [PR, PS, PL, PDV, PD] = Phase_Features(phase_sig1,  fs, Th)
%  *NOTE: While using first case the phase sequenses (phase_sig1 and 
%         phase_sig2) have to be row vectors. In secons case, the
%         phase_sig1 have to be a (2*n) matrix where the rows are phase
%         signals.
%   Single-channel : 
%   [IF, PR, PS, PL] = Phase_Features(phase_sig1,  fs, Th)
%  *NOTE: In this case the phase_sig1 is a vector containing instantaneous
%         phase signal.
% 
%             Inputs:
%                       phase_sig1 : instantaneous phase vector or matrix
%                       fs : sampling frequency (Hz)
%                       Th : Threshold value used to detect phase shift
%                            events (Radians)
%                       phase_sig2 : instantaneous phase vector
% 
%             Outputs: 
%                       PS : Phase Shift Events
%                       PL : Phase Lock Events
%                       PR : Phase Resetting Events
%                       PD : Phase Difference
%                       PDV : Phase Difference Derivatives
%                       IF : Instantaneous Frequency
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
[m1, n1] = size(phase_sig1);
if nargin==4
    phase_sig2 = varargin{1};
    [m2, n2] = size(phase_sig2);
    if m1~=1 || m2~=1
        error('***Phase Sequences Must Contain Only One Channel of Information***')
    end
    if m1~=m2 || n1~=n2
        error('***In Multi-signal Case, Phase Sequences Have to Be in Same Length***')
    end
elseif nargin>4
    error('Too Many Input Arguments: Please Check the Instructions')
end

%%-Calculating Phase features-%%
switch nargin
    case 4
        flg = 2;
        PD = phase_sig1 - phase_sig2;
        PDV = [diff(PD), 0];
        PS = zeros(m1 ,n1);
        for j=1:n1-1
            if PDV(j)>=Th || PDV(j)<=-Th
                PS(j) = 1+Th/100;
            end
        end
        PL = ~PS;
        PR = ones(m1, n1)+Th/100;    % shifting a bit to have a better view
        ind = find(PL~=0);
        for c=1:length(ind)-1
            if ind(c+1)~=ind(c)+1
                PR(1, ind(c)+1:end) = -PR(1, ind(c)+1:end);
            end
        end
        PL = PL+Th/100;
        
    case 3
        if m1==1
            flg = 1;
            PD = [diff(phase_sig1), 0];
            IF = fs*PD/2/pi;
            PDV = [diff(PD), 0];
            PS = zeros(m1 ,n1);
            for j=1:n1-1
                if PDV(j)>=Th || PDV(j)<=-Th
                    PS(j) = 1+Th/100;
                end
            end
            PL = ~PS;
            PR = ones(m1, n1)+Th/100;
            ind = find(PL~=0);
            for c=1:length(ind)-1
                if ind(c+1)~=ind(c)+1
                    PR(1, ind(c)+1:end) = -PR(1, ind(c)+1:end);
                end
            end
            PL = PL+Th/100;
        end
        
        if m1==2
            phase_sig2 = phase_sig1(2, :);
            flg = 2;
            PD = phase_sig1 - phase_sig2;
            PDV = [diff(PD), 0];
            PS = zeros(m1 ,n1);
            for j=1:n1-1
                if PDV(j)>=Th || PDV(j)<=-Th
                    PS(j) = 1+Th/100;
                end
            end
            PL = ~PS;
            PR = ones(m1, n1)+Th/100;    
            ind = find(PL~=0);
            for c=1:length(ind)-1
                if ind(c+1)~=ind(c)+1
                    PR(1, ind(c)+1:end) = -PR(1, ind(c)+1:end);
                end
            end
            PL = PL+Th/100;
            
        elseif m1>2
            error('***Input Phase Matrix Must Contain at Most TWO Channels***')
        end
        
    otherwise
        error('***Either the Fs or Threshold Value for Detecting Phase Shift Events (Th_PS) are not Specified***')
end
if flg == 2
    varargout = {PR, PS, PL, PDV, PD};
elseif flg == 1
    varargout = {IF, PR, PS, PL};
else
    error('***Something is NOT Right. Please Double-Check the Required Arguments!***')
end

end