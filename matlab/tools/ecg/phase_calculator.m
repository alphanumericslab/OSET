function [phase, phasepos] = phase_calculator(peaks)
% phase_calculator - ECG phase calculation from a given set of R-peaks.
%
% Usage:
%   [phase, phasepos] = phase_calculator(peaks)
%
% Inputs:
%   peaks: Vector of R-peak pulse train
%
% Outputs:
%   phase: The calculated phases ranging from -pi to pi. The R-peaks are located at phase = 0.
%   phasepos: The calculated phases ranging from 0 to 2*pi. The R-peaks are again located at phasepos = 0.
%
% References and Usages:
% - Sameni, R., Jutten, C., & Shamsollahi, M. B. (2008). Multichannel electrocardiogram decomposition using periodic component analysis. IEEE Transactions on Biomedical Engineering, 55(8), 1935-1940.
% - Sameni, R., Shamsollahi, M. B., Jutten, C., & Clifford, G. D. (2007). A nonlinear Bayesian filtering framework for ECG denoising. IEEE Transactions on Biomedical Engineering, 54(12), 2172-2185.
%
% Revision History:
%   2006: First release
%   2023: Renamed from deprecated version PhaseCalculation
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

phasepos = zeros(1, length(peaks));

I = find(peaks);
for i = 1:length(I)-1
    m = I(i+1) - I(i);
    phasepos(I(i)+1:I(i+1)) = 2*pi/m : 2*pi/m : 2*pi;
end
m = I(2) - I(1);
L = length(phasepos(1:I(1)));
phasepos(1:I(1)) = 2*pi - (L-1)*2*pi/m : 2*pi/m : 2*pi;

m = I(end) - I(end-1);
L = length(phasepos(I(end)+1:end));
phasepos(I(end)+1:end) = 2*pi/m : 2*pi/m : L*2*pi/m;

phasepos = mod(phasepos, 2*pi);

phase = phasepos;
I = find(phasepos > pi);
phase(I) = phasepos(I) - 2*pi;
