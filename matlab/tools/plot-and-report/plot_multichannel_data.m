function plot_multichannel_data(data, L, varargin)
% plot_multichannel_data - Plot multichannel signals in several panels and figures.
%
% Usage:
%   plot_multichannel_data(data, L, color, fs, title)
%
% Inputs:
%   data: The input data matrix
%   L: Number of panels per figure
%   color: The color of the plots: 'b', 'r', 'g', etc. (blue by default)
%   fs: Sampling rate (1 by default)
%   title: Title of the plot (empty by default)
%
% Revision History:
%   2008: First release
%   2018: Added title and sampling frequency inputs
%   2023: Renamed from deprecated version PlotECG and added sorting mode
%
% Reza Sameni, 2008-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET


% Input arguments
if nargin > 2 && ~isempty(varargin{1})
    color = varargin{1};
else
    color = 'b';
end

if nargin > 3 && ~isempty(varargin{2})
    fs = varargin{2};
    lbl = 'time(s)';
else
    fs = 1;
    lbl = 'index';
end

if nargin > 4 && ~isempty(varargin{3})
    ttl = varargin{3};
else
    ttl = '';
end

if size(data, 1) > size(data, 2)
    data = data';
end

L1 = size(data, 1);
L2 = size(data, 2);

t = (0:L2-1) / fs;

for i = 1:L1
    if mod(i, L) == 1 || L == 1
        figure;
    end
    subplot(L, 1, mod(i-1, L)+1);
    plot(t, data(i, :), color);
    ylabel(num2str(i));
    grid;
    if mod(i, L) == 1 || L == 1
        title(ttl);
    end
    if mod(i, L) == 0 || L == 1
        xlabel(lbl);
    end
end
