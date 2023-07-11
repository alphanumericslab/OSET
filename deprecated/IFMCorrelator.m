function f = IFMCorrelator(x, f0, BW, tau)
% IFMCorrelator has been deprecated. Use ifm_correlator instead.
warning('IFMCorrelator has been deprecated. Use ifm_correlator instead.');
f = ifm_correlator(x, f0, BW, tau);