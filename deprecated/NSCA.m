function [y, W, A, B, C] = NSCA(x, I, J)
% NSCA has been deprecated. Use nonstationary_component_analysis instead.
warning('NSCA has been deprecated. Use nonstationary_component_analysis instead.');
[y, W, A, B, C] = nonstationary_component_analysis(x, I, J);