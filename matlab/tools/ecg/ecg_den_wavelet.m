function [y_med, y_mean] = ecg_den_wavelet(x, varargin)
% ecg_den_wavelet - A wavelet-based ECG denoiser using majority voting or
%   averaging between multiple wavelet denoising configurations.
%
% Syntax:
%   [y_med, y_mean] = ecg_den_wavelet(x)
%   [y_med, y_mean] = ecg_den_wavelet(x, params)
%
% Inputs:
%   x: Matrix of input ECG data (channels x samples)
%   params: Struct containing denoising parameters (optional)
%       Parameters include:
%       - TPTR: Threshold selection rule (default: 'rigrsure')
%       - SORH: Soft or hard thresholding (default: 's')
%       - SCAL: Threshold scaling mode (default: 'sln')
%       - WLEVELS: Vector of wavelet decomposition levels (default: [5, 7, 9])
%       - WNAME: Cell array of wavelet names (default: {'coif2', 'coif3', 'coif4', 'coif5', 'sym5', 'sym6', 'sym7', 'sym8', 'db2'})
%
% Outputs:
%   y_med: Matrix of denoised ECG data using median of wavelet denoising outputs
%   y_mean: Matrix of denoised ECG data using mean of wavelet denoising outputs
%
% References:
%   The default parameters were found to outperform other wavelet
%   configurations over multiple ECG datasets. Read more:
%       R. Sameni, Online filtering using piecewise smoothness priors: Application
%       to normal and abnormal electrocardiogram denoising, Signal Processing,
%       Volume 133, 2017, Pages 52-63, https://doi.org/10.1016/j.sigpro.2016.10.019.
%
% Revision History:
%   2023: First release
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Default parameter values
TPTR = {'rigrsure'}; % {'rigrsure', 'heursure', 'sqtwolog', 'minimaxi'};
SORH = {'s'}; % {'s', 'h'};
SCAL = {'sln'}; % {'one', 'sln', 'mln'};
WLEVELS = [5, 7, 9]; % 1:10
WNAME = {'coif2', 'coif3', 'coif4', 'coif5', 'sym5', 'sym6', 'sym7', 'sym8', 'db2'}; % {'haar', 'db2', 'db3' ,'db4', 'db5', 'db6', 'db7', 'db8', 'db9', 'db10', 'db12', 'db16', 'coif1', 'coif2', 'coif3', 'coif4', 'coif5', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'bior1.5', 'bior2.6', 'bior2.8', 'bior5.5', 'bior6.8'};

% Check if denoising parameters are provided in varargin
if nargin > 1 && ~isempty(varargin{1})
    params = varargin{1};
    if isfield(params, 'TPTR')
        TPTR = params.TPTR;
    end
    if isfield(params, 'SORH')
        SORH = params.SORH;
    end
    if isfield(params, 'SCAL')
        SCAL = params.SCAL;
    end
    if isfield(params, 'WLEVELS')
        WLEVELS = params.WLEVELS;
    end
    if isfield(params, 'WNAME')
        WNAME = params.WNAME;
    end
end

y_med = zeros(size(x));
y_mean = zeros(size(x));

% Loop over all channels
for ch = 1 : size(x, 1)
    % Create a matrix to store wavelet denoising outputs for all configurations
    x_filtered_wavelet = zeros(length(TPTR), length(SORH), length(SCAL), length(WLEVELS), length(WNAME), length(x));

    % Loop through all denoising configurations
    for tptr = 1 : length(TPTR)
        for sorh = 1 : length(SORH)
            for scal = 1 : length(SCAL)
                for wlevel = 1 : length(WLEVELS)
                    for wname = 1 : length(WNAME)
                        % Apply wavelet denoising to the input signal
                        x_filtered_wavelet(tptr, sorh, scal, wlevel, wname, :) = wden(x(ch, :), TPTR{tptr}, SORH{sorh}, SCAL{scal}, WLEVELS(wlevel), WNAME{wname});
                    end
                end
            end
        end
    end

    % Stack the denoising outputs for median and mean calculations
    x_filtered_wavelet_stacked = reshape(x_filtered_wavelet, [], length(x));

    % Calculate the median and mean of the stacked denoising outputs
    y_med(ch, :) = median(x_filtered_wavelet_stacked, 1, 'omitnan');
    y_mean(ch, :) = mean(x_filtered_wavelet_stacked, 1, 'omitnan');
end
