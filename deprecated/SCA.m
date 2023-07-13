function [y,W,A] = SCA(x, fl, fu)
% SCA has been deprecated. Use spectral_component_analysis_dft instead.
warning('SCA has been deprecated. Use spectral_component_analysis_dft instead.');
[y, W, A] = spectral_component_analysis_dft(x, fl, fu);