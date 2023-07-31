function [phase_avg, freq_avg, amp_avg, analytic_sig_avg] = phase_est(sig, fs, f0, bw_base, varargin)
% 
% [phase_avg, freq_avg, amp_avg, analytic_sig_avg] = phase_est(sig, fs, f0, bw_base)
% [phase_avg, freq_avg, amp_avg, analytic_sig_avg] = phase_est(sig, fs, f0, bw_base, pertnum)
% 
% *************************************************************************
% * Instantaneous Phase Estimation By Transfer Function Perturbation (TFP)* 
% * Method (R. Sameni et. al. 2017 & E. Seraj et. al. 2017)               *
% *************************************************************************
% 
% Usage:    [phase_avg, freq_avg, amp_avg, analytic_sig_avg] = phase_est(sig, fs, f0, bw_base)
%           [phase_avg, freq_avg, amp_avg, analytic_sig_avg] = phase_est(sig, fs, f0, bw_base, pertnum)
% inputs:
%           'sig': input raw signal
%           'fs': sampling frequency
%           'f0': center frequency of the passband 
%           'bw_base': bandwidth of the frequency filter
%     (opt) 'pertnum': number of perturbations while using the TFP phase
%                      estimation method (default: 100)
% outputs:
%           'phase_avg': estimated instantaneous phase of input signal
%           'freq_avg': estimated instantaneous frequency of input signal
%           'amp_avg': estimated instantaneous envelope of input signal
%           'analytic_sig_avg': generated analytic form of input signal
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
if nargin < 4
    error('***wrong number of input arguments. Refer to Manual for details***')
elseif nargin == 4 
    pertnum = 100;
elseif nargin > 4
    if size(varargin, 2) ~= 1
        error('***an empty bracket [] Must be assigned to not-specified values***')
    else
        if isempty(varargin)
            pertnum = 100;
        else
            if iscalar(varargin)
                pertnum = varargin;
            else
                error('***number of perturbations must be a positive scalar***')
            end
        end
    end
end
if isscalar(sig)
    error('***input signal has to be a double vector or matrix***')
end
if (~(isscalar(fs) && isscalar(f0) && isscalar(bw_base)))
    error('***fs, f0 and bw (bandwidth) values have to be scalars***')
end

%% TFP initialization and parameter specification
order = 6;              % filter order
bw_base_dev = 0.1;      % bandwidth deviation
f0_dev = 1e-6;          % center frequency deviation
dither_std = 1e-4;      % dither level

%% applying the transfer function perturbation (TFP) method
analytic_sig = zeros(pertnum, length(sig));
phase = zeros(pertnum, length(sig));
freq = zeros(pertnum, length(sig));
amp = zeros(pertnum, length(sig));
for k=1:pertnum
    %%-generating the perturbations-%%
    dither_narrow_band = BPFilter5(randn(1, length(sig)), f0/fs, bw_base/fs, order);
    bw = bw_base + bw_base_dev*rand;
    f = f0 + f0_dev*(2*rand - 1);

    %%-analytic form representation-%%
    analytic_sig(k, :) = hilbert( BPFilter5(sig , f/fs, bw/fs, order) +...
        dither_std*dither_narrow_band/std(dither_narrow_band));

    %%-instantaneous parameter estimation-%%
    phase(k, :) = unwrap(atan2(imag(analytic_sig(k, :)), real(analytic_sig(k, :))));
    freq(k, :) = fs*diff([phase(k, 1) phase(k, :)])/(2*pi);
    amp(k, :) = abs(real(analytic_sig(k, :)) + 1i*imag(analytic_sig(k, :)));
end
analytic_sig_avg = mean(analytic_sig);
phase_avg = mean(phase);
freq_avg = mean(freq);
amp_avg = mean(amp);

end