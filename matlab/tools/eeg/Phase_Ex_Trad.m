function [phase, inst_freq, amp] = Phase_Ex_Trad(x1, fs, varargin)
%
% [phase, inst_freq, amp] = Phase_Ex_Trad(x1, fs)
% [phase, inst_freq, amp] = Phase_Ex_Trad(x1, fs, WS, max_f, ndft)
% 
% *************************************************************************
% * Instantaneous Phase (IP) estimation using the conventional analytic   *
% * representation approach through FIR filtering and Hilbert Transform.  *
% * Refer to the User Guide for further details.                          *
% *************************************************************************
% 
%             Inputs:
%                       x1 : input signal
%                       fs : sampling frequency (Hz)
%            (optional) WS : Stop-band frequency (Hz) >> (BW = WS/2)
%            (optional) max_f : Maximum frequency to extract its
%                               corresponding phase information (the IPs
%                               will be estimated for frequency components
%                               from 1(Hz) to max_f(Hz)) *** Refer to
%                               User Guide for more details.
%            (optional) ndft: number of frequency bins
% 
%             Defaults:
%                       WS = 1(Hz)
%                       max_f = 30(Hz)
%                       ndft = 100
% 
%             *NOTE: While specifying a value to one of the parameters
%                    having default values, an empty bracket [] must be
%                    used for non-specified parameters. If you're not using
%                    these parameters, empty bracket is not required.
%             Outputs: 
%                       phase : Instantaneous Phase matrix for frequency
%                               components from 1(Hz) to max_f(Hz) ***
%                               Respectively from first to last ROWS.
%                       inst_freq : Instantaneous Frequency matrix for 
%                                   frequency components from 1(Hz) to 
%                                   max_f(Hz) *** Respectively from first 
%                                   to last ROWS.
%                       amp : Instantaneous Amplitude matrix for frequency
%                             components from 1(Hz) to max_f(Hz) ***
%                             Respectively from first to last ROWS. 
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
x1 = x1(:)';
if nargin<2
    error('***At Least the Sampling Frequency (Fs) Have to Be Specified***')
elseif nargin==2
    f = linspace(0, 1, 6);
    max_f = 30;
    ndft = 100;
elseif nargin>2
    if nargin~=5
        error('***An Empty Bracket [] Must be Used for the Values NOT Being Specified***')
    end
    if isempty(varargin{1})
        f = linspace(0, 1, 6);
    else
        f = linspace(0, varargin{1}, 6);
    end
    if isempty(varargin{2})
        max_f = 30;
    else
        max_f = varargin{2};
    end
    if isempty(varargin{3})
        ndft = 100;
    else
        ndft = varargin{3};
    end
end

%%-Parks-McClellan optimal FIR filter Parameters and design-%%
tic
rp = 0.1;
rs = 60;
a = [0 1 1 0];
dev = [0.01 (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20) 0.01];
[n, fo, ao, w] = firpmord(f, a, dev, fs);
b = firpm(n, fo, ao, w);
% fvtool(b, 1, 512, fs)        % uncomment to visualize the designed filter
disp('*******************************************************************')
t1 = toc;
disp('>>>Elapsed Time for Designing the Filter in Traditional Method is:')
disp(t1)

%%-Frequency filtering and IP estimation-%%
x1_analytic = zeros(ndft, length(x1));
x1_analytic_phase = zeros(ndft, length(x1));
for i = 1:ndft,
    ff = (fs/2)*(i - 1)/ndft;
    xx1 = exp(-1j*(2*pi*ff/fs*(0:length(x1)-1))).*x1;
    
    yy1 = filter(b, 1, xx1);
    
    x1_analytic(i, :) = exp(1j*(2*pi*ff/fs*(0:length(x1)-1))).*yy1;
    
    x1_analytic(i, :) = (x1_analytic(i, :) + conj(x1_analytic(i, :)))/2;
    
    x1_analytic_phase(i, :) = hilbert(x1_analytic(i, :));
end
amp = zeros(max_f, length(x1));
phase = zeros(max_f, length(x1));
phase1 = atan2(imag(x1_analytic_phase), real(x1_analytic_phase));
phase1 = unwrap(phase1, [], 2);
for ii=1:max_f
    ch = round(ndft*ii/(fs/2)) + 1;
    amp(ii, :) = abs(real(x1_analytic_phase(ch, :)) + 1i*imag(x1_analytic_phase(ch, :)));
    phase(ii, :) = phase1(ch, :);
end
inst_freq = diff(phase, [], 2);
inst_freq = [inst_freq, zeros(max_f, 1)];
inst_freq = fs*inst_freq/2/pi;

t2 = toc;
disp('>>>Elapsed Time to Fulfill the Procedure in Traditional Method is:')
disp(t2-t1)
disp('*******************************************************************')

end