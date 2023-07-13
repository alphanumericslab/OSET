function [y, W, A] = SCAmodified(x, fl, fu)
% SCAmodified has been deprecated. Use spectral_component_analysis_bp instead.
warning('SCAmodified has been deprecated. Use spectral_component_analysis_bp instead.');
[y, W, A] = spectral_component_analysis_bp(x, fl, fu, 'DFT_FILTER');