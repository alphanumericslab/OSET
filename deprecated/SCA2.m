function [y, W, A] = SCA2(x, f0, bw, order)
% SCA2 has been deprecated. Use spectral_component_analysis_ma instead.
warning('SCA2 has been deprecated. Use spectral_component_analysis_ma instead.');
[y, W, A] = spectral_component_analysis_ma(x, f0, bw, order);