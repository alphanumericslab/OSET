function [y, W, A] = SCA_NSCA(x, f1, f2, I, order)
% SCA_NSCA has been deprecated. Use sca_nsca_component_analysis instead.
warning('SCA_NSCA has been deprecated. Use sca_nsca_component_analysis instead.');
[y, W, A] = sca_nsca_component_analysis(x, f1, f2, I, order);