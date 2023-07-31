% *************************************************************************
% * This code implements and illustrates three distributions, namely      *
% * 'Rayleigh', 'Rician' and Phase distributions for envelope and phase   *
% * introduced in [2].                                                    *
% * >>> Refer to the user manual and reference [2] for more detaiils.     *
% *************************************************************************
% 
% Dependencies: -The Cerebral Signal Phase Analysis Toolbox of Open Source 
%                Electrophysiological Toolbox
% 
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

close all
clear
clc

%%-initializing the parameters
nvar = [0.2 1 2.5];
nvar2 = [1 2 5 10 20 100];
X = 4.5;
r = 0:.001:10;
A = r;
phi = -pi:.001:pi;
fntname = 'Helvetica';          %'Times New Roman';

%%-Rayleigh distribution
RayleighD = zeros(length(nvar), length(r));
for k = 1:length(nvar),
    RayleighD(k, :) = r./nvar(k).*exp(-r.^2./(2*nvar(k)));
    AreaRayleighD = sum(RayleighD(k, :))*(r(2) - r(1))    % checking the area to not exceed "1"
end

%%-Rician distribution
RicianD = zeros(length(nvar), length(A));
for k = 1:length(nvar),
    a = X/sqrt(2*nvar(k));
    snr = 10*log10(a^2);
    RicianD(k, :) = A./nvar(k).*exp(-(A.^2 + X^2)./(2*nvar(k))).*besseli(0, X*A/nvar(k));
    AreaRicianD = sum(RicianD(k, :))*(A(2) - A(1))    % checking the area to not exceed "1"
end

%%-Phase distribution
PhaseD = zeros(length(nvar2), length(phi));
for k = 1:length(nvar2),
    a = X/sqrt(2*nvar2(k));
    snr = 10*log10(a^2);
    PhaseD(k, :) = (exp(-a^2) + a*sqrt(pi)*cos(phi).*exp(-a^2*sin(phi).^2).*erfc(-a*cos(phi)))/2/pi;
    AreaPhaseD = sum(PhaseD(k, :))*(phi(2) - phi(1))    % checking the area to not exceed "1"
end

%%-visualizing the distributions
[maxval, maxi] = max(RayleighD, [], 2);
figure
plot(r, RayleighD', 'linewidth', 2);
set(gca, 'fontsize', 16, 'fontname', fntname);
text(r(maxi), maxval+0.05, {'0.2', '1', '2.5'}, 'fontsize', 16, 'fontname', fntname);
xlabel('r_n', 'fontname', fntname);
ylabel('f(r_n|\sigma_n^2) and f(A_n|X_n,\sigma_n^2)', 'fontname', fntname);
grid
hold on
[maxval, maxi] = max(RicianD, [], 2);
plot(r, RicianD', 'linewidth', 2);
text(r(maxi), maxval+0.05, {'0.2', '1', '2.5'}, 'fontsize', 16, 'fontname', fntname);
axis([min(r) max(r) 0 1.5]);
text([1 5], [1.3 0.8], {'Background EEG', 'Foreground + Background'}, 'fontsize', 14, 'fontname', fntname);
text([1 5], [1.2 0.7], {'(Rayleigh)', 'EEG (Rician)'}, 'fontsize', 14, 'fontname', fntname);
line([3 3], [0 0.5], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
text(2.8, 0.55, 'th', 'fontsize', 14, 'fontname', fntname);

figure
plot(phi, PhaseD', 'linewidth', 2);
set(gca, 'fontsize', 16);
xlabel('\Delta\phi_n', 'fontname', fntname);
ylabel('f(\Delta\phi_n|X_n,\sigma_n^2)', 'fontname', fntname);
grid
a = axis;
a(1) = -3;
a(2) = 3;
axis(a);
