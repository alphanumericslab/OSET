function [phase_avg, freq_avg, amp_avg] = Phase_Ex_TFP(x1, fs, varargin)
% 
% [phase_avg, freq_avg, amp_avg] = Phase_Ex_TFP(x1, fs)
% [phase_avg, freq_avg, amp_avg] = Phase_Ex_TFP(x1, fs, f0, bw_base, f0_dev, bw_base_dev, dither_std, pertnum)
% 
% *************************************************************************
% * Instantaneous Phase estimation using the Transfer-Function            *
% * Perturbation Phase estimation method (TFP) [2] through analytic       *
% * representation, IIR filters and forward-backward filtering.           * 
% * >>> Refer to the User Guide for further details.                      *
% *************************************************************************
% 
%             Inputs:
%                       x1: input signal
%                       fs: sampling frequency (Hz)
%            (optional) f0: center frequency of the frequency filter (Hz)
%            (optional) bw_base: bandwidth of the frequency filter (Hz)
%            (optional) f0_dev: center frequency deviation range (Hz)
%            (optional) bw_base_dev: bandwidth deviation range (Hz)
%            (optional) dither_std: dither noise level
%            (optional) pertnum: number of attempts for perturbing filter's
%                                transfer function
% 
%             Defaults:
%                       f0 = 10(Hz)
%                       bw_base = 4(Hz) *** alpha rhythms of EEG
%                       f0_dev = 1e-6 (Hz)
%                       bw_base_dev = 0.1 (Hz)
%                       dither_std = 1e-4
%                       pertnum = 100
% 
%             *NOTE: While specifying a value to one of the parameters
%                    having default values, an empty bracket [] must be
%                    used for non-specified parameters. If you're not using
%                    these parameters, empty bracket is not required.
%             Outputs: 
%                       phase_avg: Instantaneous Phase matrix
%                       freq_avg: Instantaneous Frequency matrix
%                       amp_avg: Instantaneous Amplitude matrix 
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
    f0 = 10;
    bw_base = 4;
    f0_dev = 1e-6;
    bw_base_dev = 0.1;
    dither_std = 1e-4;
    pertnum = 100;
elseif nargin>2
    if nargin~=8
        error('***Wrong Number of Input Arguments: An Empty Bracket [] Must be Used for the Values NOT Being Specified***')
    end
    if isempty(varargin{1})
        f0 = 10;
    else
        f0 = varargin{1};
    end
    if isempty(varargin{2})
        bw_base = 4;
    else
        bw_base = varargin{2};
    end
    if isempty(varargin{3})
        f0_dev = 1e-6;
    else
        f0_dev = varargin{3};
    end
    if isempty(varargin{4})
        bw_base_dev = 0.1;
    else
        bw_base_dev = varargin{4};
    end
    if isempty(varargin{5})
        dither_std = 1e-4;
    else
        dither_std = varargin{5};
    end
    if isempty(varargin{6})
        pertnum = 100;
    else
        pertnum = varargin{6};
    end
end

% applying the transfer function perturbation (TFP) method
order = 6;              % filter order
analytic_sig = zeros(pertnum, length(x1));
phase = zeros(pertnum, length(x1));
freq = zeros(pertnum, length(x1));
amp = zeros(pertnum, length(x1));
bw = zeros(1, pertnum);
f = zeros(1, pertnum);
for k=1:pertnum
    %%-generating the perturbations-%%
    dither_narrow_band = BPFilter5(randn(1, length(x1)), f0/fs, bw_base/fs, order);
    bw(k) = bw_base + bw_base_dev*rand;
    f(k) = f0 + f0_dev*(2*rand - 1);

    %%-analytic form representation-%%
    analytic_sig(k, :) = hilbert( BPFilter5(x1 , f(k)/fs, bw(k)/fs, order) +...
        dither_std*dither_narrow_band/std(dither_narrow_band));

    %%-instantaneous parameter estimation-%%
    phase(k, :) = unwrap(atan2(imag(analytic_sig(k, :)), real(analytic_sig(k, :))));
    freq(k, :) = fs*diff([phase(k, 1) phase(k, :)])/(2*pi);
    amp(k, :) = abs(real(analytic_sig(k, :)) + 1i*imag(analytic_sig(k, :)));
end
% analytic_sig_avg = mean(analytic_sig);
phase_avg = mean(phase);
freq_avg = mean(freq);
amp_avg = mean(amp);

end