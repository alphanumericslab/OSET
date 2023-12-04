function phase = phase_shifter(phasein, teta)
% phase_shifter - Phase shifter; used for shifting the ECG phase signal
%
% Usage:
%   phase = phase_shifter(phasein, teta)
%
% Inputs:
%   phasein: Calculated ECG phase
%   teta: Desired phase shift. teta>0 and teta<0 corresponds with phase leads and phase lags, respectively.
%
% Output:
%   phase: The shifted phase.
% 
% References and usages:
% - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2008). Multichannel electrocardiogram decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering, 55(8), 1935-1940.
% - Sameni, R., Shamsollahi, M. B., Jutten, C., & Clifford, G. D. (2007). A nonlinear Bayesian filtering framework for ECG denoising. IEEE Transactions on Biomedical Engineering, 54(12), 2172-2185.
%
% Revision History:
%   2006: First release
%   2023: Renamed from deprecated version PhaseShifting
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

phase = phasein + teta;
phase = mod(phase + pi, 2 * pi) - pi;
