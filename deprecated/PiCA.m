function [y, W, A] = PiCA(x, peaks1, varargin)
% PiCA has been deprecated. Use periodic_component_analysis instead.
warning('PiCA has been deprecated. Use periodic_component_analysis instead.');
[y, W, A] = periodic_component_analysis(x, peaks1, varargin{:});
