function PlotECG(data, L, varargin)
% PlotECG has been deprecated. Use plot_multichannel_data instead.
warning('PlotECG has been deprecated. Use plot_multichannel_data instead.');
plot_multichannel_data(data, L, varargin{:});