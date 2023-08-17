function [dipole, teta] = ecg_dipole_gen_direct_sum(N, fs, f, alphai, bi, tetai, teta0)
% 
% ecg_dipole_gen_direct_sum - Synthetic cardiac dipole generator using the direct implementation
%   of the sum of Gaussian model. Performs identical to thedipole_gen_state_space.
%   Refer to references for further details.
%
% Usage:
%   [dipole, teta] = ecg_dipole_gen_direct_sum(N, fs, f, alphai, bi, tetai, teta0);
% 
% Inputs:
%   N: Signal length (in samples)
%   fs: Sampling rate
%   f: Heart rate (Hz)
%   alphai: Structure containing the amplitudes of Gaussian functions used
%       for modeling the x, y, and z coordinates of the cardiac dipole
%   bi: Structure containing the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   tetai: Structure containing the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac dipole
%   teta0: Initial phase of the synthetic dipole
%
% Outputs:
%   dipole: Structure containing the x, y, and z coordinates of the cardiac dipole
%   teta: Vector containing the dipole phase
%
% References:
%   - Sameni, R., Clifford, G. D., Jutten, C., & Shamsollahi, M. B. (2007).
%     Multichannel ECG and Noise Modeling: Application to Maternal and Fetal
%     ECG Signals. In EURASIP Journal on Advances in Signal Processing (Vol.
%     2007, Issue 1). Springer Science and Business Media LLC.
%     https://doi.org/10.1155/2007/43407
%   - McSharry, P. E., Clifford, G. D., Tarassenko, L., & Smith, L. A.
%     (2003). A dynamical model for generating synthetic electrocardiogram
%     signals. In IEEE Transactions on Biomedical Engineering (Vol. 50, Issue
%     3, pp. 289–294). Institute of Electrical and Electronics Engineers
%     (IEEE). https://doi.org/10.1109/tbme.2003.808805
% 
% Revision History:
%   2006: First release
%   2023: Renamed from the deprecated version DipoleGenerator2

% Constants
w = 2 * pi * f;
dt = 1 / fs;

% Compute dipole phase using cumulative sum
teta = cumsum([teta0 + pi, w * dt * ones(1, N - 1)]);
teta = mod(teta, 2 * pi) - pi;

% Calculate phase differences
dtetaix = mod(teta(ones(length(tetai.x), 1), :)' - tetai.x(ones(1, N), :) + pi, 2 * pi) - pi;
dtetaiy = mod(teta(ones(length(tetai.y), 1), :)' - tetai.y(ones(1, N), :) + pi, 2 * pi) - pi;
dtetaiz = mod(teta(ones(length(tetai.z), 1), :)' - tetai.z(ones(1, N), :) + pi, 2 * pi) - pi;

% Compute x, y, and z coordinates of the dipole
X = sum(alphai.x(ones(1, N), :) .* exp(-dtetaix .^ 2 ./ (2 * bi.x(ones(1, N), :) .^ 2)), 2);
Y = sum(alphai.y(ones(1, N), :) .* exp(-dtetaiy .^ 2 ./ (2 * bi.y(ones(1, N), :) .^ 2)), 2);
Z = sum(alphai.z(ones(1, N), :) .* exp(-dtetaiz .^ 2 ./ (2 * bi.z(ones(1, N), :) .^ 2)), 2);

% Store x, y, and z coordinates in the dipole structure
dipole.x = X';
dipole.y = Y';
dipole.z = Z';
