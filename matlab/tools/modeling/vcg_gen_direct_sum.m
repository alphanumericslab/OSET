function [vcg, teta] = vcg_gen_direct_sum(N, fs, f, alpha, b, theta, theta0)
% 
% vcg_gen_direct_sum - Synthetic cardiac dipole/vectorcardiogram (VCG) generator using the direct implementation
%   of the sum of Gaussian model. Performs identical to thedipole_gen_state_space.
%   Refer to references for further details.
%
% Usage:
%   [vcg, teta] = vcg_gen_direct_sum(N, fs, f, alpha, b, theta, theta0);
% 
% Inputs:
%   N: Signal length (in samples)
%   fs: Sampling rate
%   f: Heart rate (Hz)
%   alpha: Structure containing the amplitudes of Gaussian functions used
%       for modeling the x, y, and z coordinates of the cardiac VCG/dipole
%   b: Structure containing the widths of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac VCG/dipole
%   theta: Structure containing the phase of Gaussian functions used for
%       modeling the x, y, and z coordinates of the cardiac VCG/dipole
%   theta0: Initial phase of the synthetic VCG/dipole
%
% Outputs:
%   vcg: Structure containing the x, y, and z coordinates of the cardiac VCG/dipole
%   teta: Vector containing the VCG phase
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

% Compute VCG phase using cumulative sum
teta = cumsum([theta0 + pi, w * dt * ones(1, N - 1)]);
teta = mod(teta, 2 * pi) - pi;

% Calculate phase differences
dtetaix = mod(teta(ones(length(theta.x), 1), :)' - theta.x(ones(1, N), :) + pi, 2 * pi) - pi;
dtetaiy = mod(teta(ones(length(theta.y), 1), :)' - theta.y(ones(1, N), :) + pi, 2 * pi) - pi;
dtetaiz = mod(teta(ones(length(theta.z), 1), :)' - theta.z(ones(1, N), :) + pi, 2 * pi) - pi;

% Compute x, y, and z coordinates of the VCG
X = sum(alpha.x(ones(1, N), :) .* exp(-dtetaix .^ 2 ./ (2 * b.x(ones(1, N), :) .^ 2)), 2);
Y = sum(alpha.y(ones(1, N), :) .* exp(-dtetaiy .^ 2 ./ (2 * b.y(ones(1, N), :) .^ 2)), 2);
Z = sum(alpha.z(ones(1, N), :) .* exp(-dtetaiz .^ 2 ./ (2 * b.z(ones(1, N), :) .^ 2)), 2);

% Store x, y, and z coordinates in the VCG structure
vcg.x = X';
vcg.y = Y';
vcg.z = Z';
