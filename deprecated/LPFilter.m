function y = LPFilter(x, fc)
% Deprecated: LPFilter is deprecated. Use lp_filter_zero_phase instead.
warning('Deprecated: LPFilter is deprecated. Use lp_filter_zero_phase instead.');
y = lp_filter_zero_phase(x, fc);