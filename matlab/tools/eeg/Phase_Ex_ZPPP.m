function [phase, inst_freq, amp] = Phase_Ex_ZPPP(x1, fs, varargin)
% 
% [phase, inst_freq, amp] = Phase_Ex_ZPPP(x1, fs)
% [phase, inst_freq, amp] = Phase_Ex_ZPPP(x1, fs, WS, max_f, ndft, pertnum)
% 
% *************************************************************************
% * Instantaneous Phase (IP) estimation using the Zer-Pole Perturbation   *
% * Phase estimation method (Z-PPP) [1] through analytic representation,  *
% * IIR filters and forward-backward filtering.                           *
% * >>> Refer to the User Guide for further details.                      *
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
%            (optional) pertnum: number of attempts for perturbing filter's
%                                zeros and poles
% 
%             Defaults:
%                       WS = 1(Hz)
%                       max_f = 30(Hz)
%                       ndft = 100
%                       pertnum = 100
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
    f = linspace(1/2, 1, 2);
    max_f = 30;
    ndft = 100;
    pertnum = 100;
elseif nargin>2
    if nargin~=6
        error('***An Empty Bracket [] Must be Used for the Values NOT Being Specified***')
    end
    if isempty(varargin{1})
        f = linspace(1/2, 1, 2);
    else
        f = linspace(varargin{1}/2, varargin{1}, 2);
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
    if isempty(varargin{4})
        pertnum = 100;
    else
        pertnum = varargin{4};
    end
end

%%-IIR filter parameters and design-%%
tic
rp = 0.1;
rs = 70;
wp = f(1)/(fs/2);
ws = f(2)/(fs/2);
[iir_ord, cutoff] = ellipord(wp, ws, rp, rs);
[b, a] = ellip(iir_ord, rp, rs, cutoff);
% fvtool(b, a, 512, fs)        % uncomment to visualize the designed filter
disp('*******************************************************************')
t1 = toc;
disp('>>>Elapsed Time for Designing the Filter in ZPPP Method is:')
disp(t1)

%%-Perturbing Zeros & Poles-%%
[b, a] = eqtflength(b, a);
[z, p, k] = tf2zp(b, a);
z = conj(z');
p = conj(p');

kk = 1;
filt_poles = [p; zeros(pertnum, length(p))];
filt_zeros = [z; zeros(pertnum, length(z))];
real_p = zeros(1, length(p));
imag_p = zeros(1, length(p));
new_p = zeros(1, length(p));
theta_new = zeros(1, length(z));
while (kk <= pertnum)
    for ii=1:length(p)
        if imag(p(ii))>0
            real_p(:, ii) = real(p(ii)) - rand/1e3;
            imag_p(:, ii) = imag(p(ii)) - rand/8e2;
            new_p(:, ii) = real_p(:, ii) + 1i*imag_p(:, ii);
        else
            new_p(:, ii) = conj(new_p(:, ii-1));
        end
    end

    [theta_z, rho_z] = cart2pol(real(z), imag(z));
    for jj=1:length(z)
        if theta_z(jj)>0
            theta_new(:, jj) = theta_z(jj) - rand/3e2;
        else
            theta_new(:, jj) = -theta_new(:, jj-1);
        end
    end
    [x_z, y_z] = pol2cart(theta_new, rho_z);
    new_z = x_z + 1i*y_z;
   
    kk = kk + 1;
    filt_poles(kk, :) = new_p;
    filt_zeros(kk, :) = new_z;
end

%%-Frequency filtering and phase calculation-%%
x1_analytic = zeros(ndft, length(x1));
x1_analytic_phase = zeros(ndft, length(x1));
x1_analytic_phase_pertrbd = cell(1, pertnum);
for j = 1:pertnum+1
    for i = 1:ndft,
        ff = (fs/2)*(i - 1)/ndft;
        xx1 = exp(-1j*(2*pi*ff/fs*(0:length(x1)-1))).*x1;
        
        [num, denum] = zp2tf(conj(filt_zeros(j, :)'), conj(filt_poles(j, :)'), k);
        
        yy1 = filtfilt(num, denum, xx1);

        x1_analytic(i, :) = exp(1j*(2*pi*ff/fs*(0:length(x1)-1))).*yy1;

        x1_analytic_phase(i, :) = x1_analytic(i, :);
    end
    x1_analytic_phase_pertrbd{j} = x1_analytic_phase;
end
amp_pertrbd = cell(1, pertnum+1);
phase1_pertrbd = cell(1, pertnum+1);
amp_avg = cell(1, max_f);
phase_avg = cell(1, max_f);
for cc=1:pertnum+1
    x1_phase = cell2mat(x1_analytic_phase_pertrbd(cc));
    phase1 = atan2(imag(x1_phase), real(x1_phase));
    phase1 = unwrap(phase1, [], 2);
    for ll=1:max_f
        ch = round(ndft*ll/(fs/2)) + 1;
        amp_pertrbd{cc}(ll, :) = abs(real(x1_phase(ch, :)) + 1i*imag(x1_phase(ch, :)));
        phase1_pertrbd{cc}(ll, :) = phase1(ch, :);
        amp_avg{ll}(cc, :) = amp_pertrbd{cc}(ll, :);
        phase_avg{ll}(cc, :) = phase1_pertrbd{cc}(ll, :);
    end
end
amp = zeros(max_f, length(x1));
phase = zeros(max_f, length(x1));
for rr=1:max_f
    amp(rr, :) = mean(amp_avg{rr});
    phase(rr, :) = mean(phase_avg{rr});
end
inst_freq = diff(phase, [], 2);
inst_freq = [inst_freq, zeros(max_f, 1)];
inst_freq = fs*inst_freq/2/pi;

t2 = toc;
disp('>>>Elapsed Time to Fulfill the Procedure in ZPPP Method is:')
disp(t2-t1)
disp('*******************************************************************')

end