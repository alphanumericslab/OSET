function y = LPFilter(x, fc)
% LPFilter has been deprecated. Use lp_filter_zero_phase instead.
warning('LPFilter has been deprecated. Use lp_filter_zero_phase instead.');
y = lp_filter_zero_phase(x, fc);